#include "GanonBuild.hpp"

#include <robin_hood.h>

#include <defaults/defaults.hpp>
#include <utils/IBFConfig.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <cinttypes>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace GanonBuild
{

namespace detail
{


// 定义未压缩的数据布局的交错布隆过滤器类型
typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > TIBF;

// 定义哈希计数的无序映射类型
typedef robin_hood::unordered_map< std::string, uint64_t > THashesCount;

// 定义bin映射哈希的无序映射类型
// 这里的映射从uint64_t类型的哈希值到一个包含目标字符串和两个uint64_t类型整数的元组
typedef robin_hood::unordered_map< uint64_t, std::tuple< std::string, uint64_t, uint64_t > > TBinMapHash;

// 定义一个结构体来表示输入文件的映射
struct InputFileMap
{
    std::string                target;  // 目标名称
    std::vector< std::string > files;   // 文件列表
};



// 定义一个结构体来表示统计的总量
struct Total
{
    uint64_t files             = 0; // 文件数
    uint64_t invalid_files     = 0; // 无效文件数
    uint64_t sequences         = 0; // 序列数
    uint64_t skipped_sequences = 0; // 跳过的序列数
    uint64_t length_bp         = 0; // 序列长度（以碱基对为单位）
};

// 定义一个结构体来表示统计信息
struct Stats
{
    Total total;    // 总量
    void  add_totals( std::vector< Total > const& totals )
    {
        // 将多个线程的总量添加到统计信息中
        for ( auto const& t : totals )
        {
            total.files += t.files;     // 累加文件数
            total.invalid_files += t.invalid_files;     // 累加无效文件数
            total.sequences += t.sequences;     // 累加序列数
            total.skipped_sequences += t.skipped_sequences;     // 累加跳过的序列数
            total.length_bp += t.length_bp;     // 累加序列长度
        }
    }
};

inline std::string get_seqid( std::string header )
{
    /*
     * 返回 FASTA 标头中的序列 ID（直到第一个空格之前的所有内容）
     */
    return header.substr( 0, header.find( ' ' ) );
}

robin_hood::unordered_map< std::string, std::vector< std::string > > parse_input_file( const std::string& input_file,
                                                                                       THashesCount&      hashes_count,
                                                                                       bool               quiet,
                                                                                       Stats&             stats )
{
    /*
     * 解析输入文件的函数 -> 一个包含字段的表格文件: 文件 [<tab> 目标 <tab> 序列 ID]
     * 返回一个映射 {目标: {文件: [序列 ID]}}
     * 如果提供了第三行的序列 ID，则 seqids 是一组序列 ID
     * 如果是文件解析，则每个文件只有一个目标，seqids 为空
     */
    robin_hood::unordered_map< std::string, std::vector< std::string > > input_map;
    std::string                                                          line;
    robin_hood::unordered_set< std::string >                             files;
    std::ifstream                                                        infile( input_file );
    // 逐行读取输入文件
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        // 使用制表符分割每行字段
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // file [<tab> 目标 <tab> 序列 ID]

        const std::string file = fields[0];
        files.insert( file );
        // 检查文件是否存在或文件是否为空
        if ( !std::filesystem::exists( file ) || std::filesystem::file_size( file ) == 0 )
        {
            if ( !quiet )
                std::cerr << "WARNING: input file not found/empty: " << file << std::endl;
            stats.total.invalid_files++;
            continue;
        }

        if ( fields.size() == 1 )
        {
            // target 是文件本身（仅文件名，没有路径）
            const auto target = std::filesystem::path( file ).filename();
            input_map[target].push_back( file );
            hashes_count[target] = 0;
        }
        else if ( fields.size() == 2 )
        {
            // 文件中提供了 target
            input_map[fields[1]].push_back( file );
            hashes_count[fields[1]] = 0;
        }
    }
    stats.total.files = files.size();

    return input_map;
}


void store_hashes( const std::string                            target,
                   const robin_hood::unordered_set< uint64_t >& hashes,
                   const std::string                            tmp_output_folder )
{
    /*
     * 将集合中的哈希值存储到指定文件夹中的磁盘上
     * 如果未指定文件夹，则存储在当前文件夹(".")中
     */
    // 创建输出文件路径
    std::filesystem::path outf{ tmp_output_folder };
    outf += target + ".min";    // 生成文件名为 "target.min"
    // 以二进制和追加模式打开输出文件
    std::ofstream outfile{ outf, std::ios::binary | std::ios::app };
    // 将每个哈希值写入文件
    for ( auto&& h : hashes )
    {
        outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
    }
    // 关闭文件
    outfile.close();
}


std::vector< uint64_t > load_hashes( std::string file )
{
    /*
     * 从磁盘加载哈希文件并返回一个包含这些哈希值的向量
     */
    uint64_t                hash;
    std::vector< uint64_t > hashes;
    std::ifstream           infile{ file, std::ios::binary };
    // 读取文件中的所有哈希值并存入向量中
    while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
        hashes.push_back( hash );
    return hashes;
}

void delete_hashes( const THashesCount& hashes_count, const std::string tmp_output_folder )
{
    /*
     * 从磁盘删除哈希文件
     */
    for ( auto const& [target, cnt] : hashes_count )
    {
        std::filesystem::path outf{ tmp_output_folder };
        outf += target + ".min";
        if ( std::filesystem::exists( outf ) )
            std::filesystem::remove( outf );
    }
}

void count_hashes( SafeQueue< InputFileMap >& ifm_queue,
                   THashesCount&              hashes_count,
                   const GanonBuild::Config&  config,
                   Total&                     total )
{

    /*
     * 该函数在多个线程中并行运行，从一个 InputFileMap 队列中读取数据
     * 计算每个目标的最小化器哈希值并将其存储到磁盘
     * 它还会保留这些计数，以便进一步确定最佳过滤器大小
     */


    // 每个线程一个视图，用于计算最小化器哈希
    auto minimiser_view = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( config.kmer_size ) } );

    while ( true )
    {
        // 等待直到有读取数据可用，或者推送结束并且队列为空
        InputFileMap ifm = ifm_queue.pop();

        // 如果在弹出后为空，退出线程
        if ( ifm.target == "" )
            break;

        robin_hood::unordered_set< uint64_t > target_hashes;

        // 针对目标的所有文件
        for ( auto& file : ifm.files )
        {
            if ( target_hashes.size() >= 2000000 )
                break;
            try
            {
                // 打开文件
                seqan3::sequence_file_input< raptor::dna4_traits,
                                             seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                    fin{ file };


                // 将文件作为目标 - 从可能包含多个序列的文件中生成所有哈希值
                // 然后进行计数和存储

                robin_hood::unordered_set< uint64_t > hashes;
                for ( auto const& [header, seq] : fin )
                {
                    if ( seq.size() < config.min_length )
                    {
                        total.skipped_sequences++;
                        continue;
                    }

                    total.sequences++;
                    total.length_bp += seq.size();
                    const auto mh = seq | minimiser_view | std::views::common;

                    for ( auto hash : mh )
                    {
                        if ( target_hashes.size() >= 2'000'000 )
                            break;
                        target_hashes.insert( hash );
                    }

                    // 达到上限后立即停止处理当前文件
                    if ( target_hashes.size() >= 2'000'000 )
                        break;

                }
                hashes_count[ifm.target] += target_hashes.size();
                detail::store_hashes( ifm.target, target_hashes, config.tmp_output_folder );
            }
            catch ( seqan3::parse_error const& e )
            {
                std::cerr << "Error parsing file [" << file << "]. " << e.what() << std::endl;
                continue;
            }
        }
    }
}

void save_filter( const GanonBuild::Config& config,
                  const TIBF&               ibf,
                  const IBFConfig&          ibf_config,
                  const THashesCount&       hashes_count,
                  const TBinMapHash&        bin_map_hash )
{
    /*
     * 将结构体和 IBF 保存到文件中。结构体被转换为标准库格式。
     * 保存内容包括：
     * 版本信息:   tuple(major, minor, patch)
     * ibf_config: 用于构建过滤器的 IBFConfig 结构体
     * hashes_count: vector<tuple(string, int)> [(目标, 计数)]
     * bin_map:     vector<tuple(int, string)> [(bin编号, 目标)]
     * ibf:         seqan3 中的交织布隆过滤器
     */
    // 创建一个二进制输出文件流
    std::ofstream               os( config.output_file, std::ios::binary );
    // 创建一个二进制输出档案
    cereal::BinaryOutputArchive archive( os );

    // 创建每个 bin 的映射并存储在标准库结构中 {bin编号: 目标}
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    for ( auto& [binno, vals] : bin_map_hash )
    {
        bin_map.push_back( { binno, std::get< 0 >( vals ) } );
    }

    // 将 hashes_count 存储在标准库结构中 {目标: 哈希计数}
    std::vector< std::tuple< std::string, uint64_t > > hashes_count_std;
    for ( auto& [target, count] : hashes_count )
    {
        hashes_count_std.push_back( { target, count } );
    }

    // 序列化保存数据
    archive( std::make_tuple( defaults::version_major, defaults::version_minor, defaults::version_patch ) );
    archive( ibf_config );
    archive( hashes_count_std );
    archive( bin_map );
    archive( ibf );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes )
{
    /*
     * 根据最大假阳性率和要插入的元素数量计算 bin 大小 (布隆过滤器)。
     */
    return std::ceil( ( n_hashes * std::log( max_fp ) ) / std::log( 1.0 / std::pow( 2, std::log( 2 ) ) ) );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes, uint8_t hash_functions )
{
    /*
     * 根据最大假阳性率、要插入的元素数量和哈希函数的数量计算 bin 大小 (布隆过滤器)。
     */
    return std::ceil( n_hashes
                      * ( -hash_functions / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) ) );
}

uint8_t hash_functions_from_ratio( uint64_t bin_size_bits, uint64_t n_hashes )
{
    /*
     * 根据 bin 大小和元素数量计算最佳哈希函数数量。
     */
    return static_cast< uint8_t >( std::log( 2 ) * ( bin_size_bits / static_cast< double >( n_hashes ) ) );
}

uint8_t get_optimal_hash_functions( uint64_t bin_size_bits,
                                    uint64_t n_hashes,
                                    uint8_t  hash_functions,
                                    uint8_t  max_hash_functions )
{
    /*
     * 该函数检查是否应该根据比率计算哈希函数数量 (hash_functions = 0)，
     * 并确保哈希函数数量在允许范围内 (1-5)。
     */

    uint8_t optimal_hash_functions = hash_functions;
    // 如果哈希函数数量为0，则根据比率计算最佳哈希函数数量
    if ( optimal_hash_functions == 0 )
        optimal_hash_functions = hash_functions_from_ratio( bin_size_bits, n_hashes );
    // 如果计算出的哈希函数数量超过最大值或为0，则将其设置为最大值
    if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
        optimal_hash_functions = max_hash_functions;

    return optimal_hash_functions;
}


uint64_t number_of_bins( THashesCount const& hashes_count, uint64_t n_hashes )
{
    /*
     * 根据元素数量计算 IBF 的 bin 数量（包括拆分的 bin）。
     */
    uint64_t n_bins = 0;
    // 遍历每个目标及其对应的哈希计数
    for ( auto const& [target, count] : hashes_count )
    {
        // 计算目标所需的 bin 数量，并累加到总 bin 数量中
        n_bins += std::ceil( count / static_cast< double >( n_hashes ) );
    }
    return n_bins;
}


double correction_rate( uint64_t max_split_bins, double max_fp, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * 计算 bin 的大小应增加的比率，以适应由将目标拆分成多个 bin 而产生的多重错误问题。
     * 基于目标在最多的 bin 中拆分的情况。
     */

    // 计算目标错误率
    double const target_fpr        = 1.0 - std::exp( std::log( 1.0 - max_fp ) / max_split_bins );
    // 计算新的 bin 大小
    size_t const new_bin_size      = bin_size( target_fpr, n_hashes, hash_functions );
    // 计算原始的 bin 大小
    size_t const original_bin_size = bin_size( max_fp, n_hashes, hash_functions );
    // 返回新的 bin 大小与原始 bin 大小的比率
    return static_cast< double >( new_bin_size ) / original_bin_size;
}


inline uint64_t optimal_bins( uint64_t n_bins )
{
    /*
     * 返回创建 IBF 的最优 bin 数量（64 的倍数）
     */
    return std::ceil( n_bins / 64.0 ) * 64;
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * 根据参数计算 bin 的理论误报率（布隆过滤器）
     */
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( n_hashes ) ) ),
                     hash_functions );
}

std::tuple< double, double > true_false_positive( THashesCount const& hashes_count,
                                                  uint64_t            max_hashes_bin,
                                                  uint64_t            bin_size_bits,
                                                  uint8_t             hash_functions )
{
    /*
     * 计算 IBF 的理论误报率（平均和最大），基于targets和split的 bins
     */
    double highest_fp = 0;
    double average_fp = 0;
    // 计算每个目标组的真实误报率，考虑拆分成多个 bins（多次测试）
    for ( auto const& [target, count] : hashes_count )
    {
        // 使用每个 bin 的平均哈希数来计算误报率
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( max_hashes_bin ) );
        // 这可能会有一个非常小的误差（四舍五入到多个 bins）
        uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

        // 当前目标的误报率，考虑拆分的 bins
        double real_fp =
            1.0 - std::pow( 1.0 - false_positive( bin_size_bits, hash_functions, n_hashes_bin ), n_bins_target );

        if ( real_fp > highest_fp )
            highest_fp = real_fp;
        average_fp += real_fp;
    }
    average_fp = average_fp / static_cast< double >( hashes_count.size() );

    return std::make_tuple( highest_fp, average_fp );
}

inline uint64_t get_max_hashes( THashesCount const& hashes_count )
{
    /*
     * 返回 `hashes_count` 中的最大哈希数
     */
    uint64_t max_hashes = 0;
    for ( auto const& [target, cnt] : hashes_count )
    {
        if ( cnt > max_hashes )
            max_hashes = cnt;
    }
    return max_hashes;
}

void optimal_hashes( double const        max_fp,
                     double const        filter_size,
                     IBFConfig&          ibf_config,
                     THashesCount const& hashes_count,
                     uint8_t const       hash_functions,
                     uint8_t const       max_hash_functions,
                     std::string const   mode )
{

    /*
     * 给定一个最大假阳性率或过滤器大小，迭代可能的 bin 容量（单个布隆过滤器），并计算相应的大小，考虑拆分的 bin
     * 选择生成"最佳"之间最小的过滤器/假阳性率和 bin 数量的参数，并填充 ibf_config 结构
     */

    // 拥有最多最小化器的目标
     uint64_t const max_hashes = get_max_hashes( hashes_count );

    // 保存平均值的最小值
    uint64_t min_filter_size = 0;
    uint64_t min_bins        = 0;
    double   min_fp          = 1;

    // 保存用于平均值的模拟值
    struct SimParam
    {
        uint64_t n_hashes;
        uint64_t n_bins;
        uint64_t filter_size_bits;
        double   fp;
    };
    std::vector< SimParam > simulations;

    // 每 100 个元素进行一次模拟
    size_t iter = 100;
    // 检查 max_hashes 是否小于迭代步长
    if ( max_hashes < iter )
        iter = max_hashes;

    // (total + 1) 用于处理零索引
    for ( size_t n = max_hashes + 1; n > iter; n -= iter )
    {
        // 要插入到一个 bin 中的元素数量
        uint64_t const n_hashes = n - 1;

        // 基于目标和元素的实际 bin 数量（不是 64 的倍数）
        uint64_t const n_bins = number_of_bins( hashes_count, n_hashes );

        int64_t bin_size_bits          = 0;
        uint8_t optimal_hash_functions = 0;
        if ( filter_size )
        {
            // 如果提供了 filter_size，简单计算
            bin_size_bits = ( filter_size / static_cast< double >( optimal_bins( n_bins ) ) ) * 8388608u;
            optimal_hash_functions =
                get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
        }
        else
        {
            if ( hash_functions == 0 )
            {
                // 首先定义大小，然后定义哈希函数数量
                bin_size_bits = bin_size( max_fp, n_hashes );
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
            }
            else
            {
                // 提供了哈希函数数量，使用它定义最佳 bin 大小
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
                bin_size_bits = bin_size( max_fp, n_hashes, optimal_hash_functions );
            }
        }

        // 目标拆分为多个 bin 的最大次数
        uint64_t const max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        // 计算 fp 或 filter_size，取决于提供的参数
        double   fp               = 0;
        uint64_t filter_size_bits = 0;
        if ( filter_size )
        {
            // 针对当前值计算假阳性率，考虑拆分的 bin
            fp =
                1 - std::pow( 1.0 - false_positive( bin_size_bits, optimal_hash_functions, n_hashes ), max_split_bins );

            // save minimal value
            // 保存最小值
            if ( fp < min_fp )
                min_fp = fp;
        }
        else
        {

            // 实际插入到每个 bin 的元素数量
            uint64_t const avg_n_hashes = std::ceil( max_hashes / static_cast< double >( max_split_bins ) );
            // 基于每个拆分 bin 的平均 n_hashes 近似真实假阳性率
            // 如果不适用，则会高估修正率，如果 bin 没有完全“满”
            double approx_fp = false_positive( bin_size_bits, optimal_hash_functions, avg_n_hashes );
            // 如果近似值较高（精度计算），则设回
            if ( approx_fp > max_fp )
                approx_fp = max_fp;

            // 基于单个目标最大拆分数量的修正率
            double const crate = correction_rate( max_split_bins, approx_fp, optimal_hash_functions, n_hashes );
            // 应用于 bin 大小
            bin_size_bits = bin_size_bits * crate;
            // 计算最终过滤器大小
            filter_size_bits = bin_size_bits * optimal_bins( n_bins );

            // 值过小或由于小的 n_hashes 或过高的 crate 导致值过大时退出循环
            if ( filter_size_bits == 0 || std::isinf( crate ) )
                break;

            // 保存最小值
            if ( filter_size_bits < min_filter_size || min_filter_size == 0 )
                min_filter_size = filter_size_bits;
        }

        // 保存模拟值
        simulations.emplace_back( SimParam{ n_hashes, n_bins, filter_size_bits, fp } );

        /*        std::cout << "n_hashes: " << n_hashes << '\t';
                std::cout << "n_bins: " << n_bins << '\t';
                std::cout << "filter_size_bits: " << filter_size_bits << '\t';
                std::cout << "fp: " << fp << '\t';
                std::cout << "hash_functions: " << unsigned( optimal_hash_functions ) << '\n';*/

        // 保存最小值
        if ( n_bins < min_bins || min_bins == 0 )
            min_bins = n_bins;
    }

    // 选择“最优”哈希作为 n_bins 和 filter_size 的调和平均数
    // 考虑它们与可能的最小值的差异
    // 如果选择了特殊模式，平均值会被一个因子偏移（对于比率来说更小更好，因此为 0.5）
    // 0 表示忽略度量并仅使用其他值（最快或最小）
    double mode_val = 1; // 1 是默认模式 avg，n_bins 和 filter_size 的调和平均数
    if ( mode == "smaller" || mode == "faster" )
        mode_val = 0.5;
    else if ( mode == "smallest" || mode == "fastest" )
        mode_val = 0;

    // 用于 fp 或 filter_size_bits，取决于提供的参数
    double var_val = 1;
    // 用于 bin 比率
    double bins_val = 1;
    if ( mode == "smaller" || mode == "smallest" )
        var_val = mode_val;
    else if ( mode == "faster" || mode == "fastest" )
        bins_val = mode_val;

    // 初始化最小平均值为0
    double min_avg = 0;
    // 遍历所有的模拟参数
    for ( auto const& params : simulations )
    {
        // 变量比率，用于比较当前参数和最小值之间的比例
        double var_ratio = 0;
        // 根据是否提供filter_size来计算变量比率
        if ( filter_size )
            // 如果提供了filter_size，则计算当前假阳性率与最小假阳性率的比率
            var_ratio = params.fp / min_fp;
        else
            // 如果没有提供filter_size，则计算当前过滤器大小与最小过滤器大小的比率
            var_ratio = params.filter_size_bits / static_cast< double >( min_filter_size );

        // 计算bins比率，当前bins数量与最小bins数量的比率
        double const bins_ratio = params.n_bins / static_cast< double >( min_bins );
        // 计算平均值，根据给定模式的值来调整权重
        // 1 + mode_val^2 确保权重不为零
        // (var_ratio * bins_ratio) 是总体比率
        // (var_val * var_ratio) + (bins_val * bins_ratio) 是加权比率
        double const avg        = ( 1 + std::pow( mode_val, 2 ) )
                           * ( ( var_ratio * bins_ratio ) / ( ( var_val * var_ratio ) + ( bins_val * bins_ratio ) ) );


        // 如果当前平均值小于最小平均值或最小平均值为0，则更新最小平均值和ibf配置
        if ( avg < min_avg || min_avg == 0 )
        {
            min_avg = avg;  // 更新最小平均值
            // 根据提供的参数更新ibf配置
            if ( filter_size )
            {
                // 如果提供了filter_size，则计算并设置ibf的bin大小和最大假阳性率
                ibf_config.bin_size_bits =
                    ( filter_size / static_cast< double >( optimal_bins( params.n_bins ) ) ) * 8388608u;
                ibf_config.max_fp = params.fp;
            }
            else
            {
                // 如果没有提供filter_size，则使用当前参数的过滤器大小和给定的最大假阳性率
                ibf_config.bin_size_bits = params.filter_size_bits / optimal_bins( params.n_bins );
                ibf_config.max_fp        = max_fp;
            }

            // 设置ibf配置的最大哈希数、bins数量和哈希函数数量
            ibf_config.max_hashes_bin = params.n_hashes;
            ibf_config.n_bins         = params.n_bins;
            ibf_config.hash_functions = get_optimal_hash_functions(
                ibf_config.bin_size_bits, params.n_hashes, hash_functions, max_hash_functions );
        }
    }
}


TBinMapHash create_bin_map_hash( IBFConfig const& ibf_config, THashesCount const& hashes_count )
{

    /*
     * 根据IBF参数和当前的哈希值创建bin编号
     * 将哈希值按平均值分配到分割的bin中
     * 输出一个包含bin编号和哈希值索引的映射，用于插入到IBF中
     */
    uint64_t    binno = 0;  // 初始化bin编号
    TBinMapHash bin_map_hash;   // 创建bin映射哈希表
    for ( auto const& [target, count] : hashes_count )
    {
        // 计算目标需要的平均bin数量
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        // 计算每个bin中包含的哈希数量
        uint64_t n_hashes_bin  = std::ceil( count / static_cast< double >( n_bins_target ) );

        // 如果计算出的哈希数量超过最大允许值，则设置为最大值
        if ( n_hashes_bin > ibf_config.max_hashes_bin )
            n_hashes_bin = ibf_config.max_hashes_bin;
        // 遍历每个目标的bin
        for ( uint64_t i = 0; i < n_bins_target; ++i )
        {
            uint64_t hashes_idx_st = i * n_hashes_bin;  // 计算当前bin的起始索引
            uint64_t hashes_idx_en = hashes_idx_st + n_hashes_bin - 1;  // 计算当前bin的结束索引
            // 如果起始索引超出总哈希数量，跳出循环
            if ( hashes_idx_st >= count )
                break;
            // 如果结束索引超出总哈希数量，设置为最大索引
            if ( hashes_idx_en >= count )
                hashes_idx_en = count - 1;
            // 将当前bin的信息插入映射哈希表中
            bin_map_hash[binno] = std::make_tuple( target, hashes_idx_st, hashes_idx_en );
            binno++;    // 增加bin编号
        }
    }
    // 断言创建的bin数量与预期相同
    assert( bin_map_hash.size() == ibf_config.n_bins );
    return bin_map_hash;    // 返回bin映射哈希表
}

void build( TIBF&                       ibf,
            std::atomic< std::size_t >& bin_batches,
            const uint64_t              max_batch,
            const uint64_t              batch_size,
            const TBinMapHash&          bin_map_hash,
            std::string                 tmp_output_folder )
{
    /*
     * 从最小化器文件和先前生成的映射构建 IBF
     */
    while ( true )
    {
        // 增加原子计数器 bin_batches 并存储当前批次编号
        uint64_t batch = bin_batches++;
        if ( batch >= max_batch )
            break;

        // 设置并检查批次的边界
        uint64_t batch_start = batch * batch_size;
        uint64_t batch_end   = batch_start + batch_size - 1;
        if ( batch_end > bin_map_hash.size() - 1 )
            batch_end = bin_map_hash.size() - 1;

        // 存储文件和哈希到一个映射中以便快速访问
        // 因为 bin 批次可能会有重复和多个文件
        robin_hood::unordered_map< std::string, std::vector< uint64_t > > target_hashes;

        // 按索引将哈希插入到 ibf 中
        for ( uint64_t binno = batch_start; binno <= batch_end; binno++ )
        {
            auto [target, hashes_idx_st, hashes_idx_en] = bin_map_hash.at( binno );
            // 只读取一次文件并存储
            if ( !target_hashes.count( target ) )
            {
                auto file             = tmp_output_folder + target + ".min";
                target_hashes[target] = load_hashes( file );
            }
            // 将哈希插入到相应的 bin 中
            for ( uint64_t pos = hashes_idx_st; pos <= hashes_idx_en; pos++ )
            {
                ibf.emplace( target_hashes[target][pos], seqan3::bin_index{ binno } );
            }
        }
    }
}

void print_stats( Stats& stats, const IBFConfig& ibf_config, const StopClock& timeGanonBuild )
{
    /*
     * 打印构建的总体统计信息
     */
    double elapsed = timeGanonBuild.elapsed();
    std::cerr << "ganon-build processed " << stats.total.sequences << " sequences / " << stats.total.files << " files ("
              << stats.total.length_bp / 1000000.0 << " Mbp) in " << elapsed << " seconds ("
              << ( stats.total.length_bp / 1000000.0 ) / ( elapsed / 60.0 ) << " Mbp/m)" << std::endl;

    // 如果存在无效文件，打印无效文件数量
    if ( stats.total.invalid_files > 0 )
        std::cerr << " - " << stats.total.invalid_files << " invalid files skipped" << std::endl;
    // 如果存在跳过的序列，打印跳过的序列数量
    if ( stats.total.skipped_sequences > 0 )
        std::cerr << " - " << stats.total.skipped_sequences << " sequences skipped" << std::endl;

    // 打印最大和平均的假阳性率，精确到4位小数
    std::cerr << std::fixed << std::setprecision( 4 ) << " - max. false positive: " << ibf_config.true_max_fp;
    std::cerr << std::fixed << std::setprecision( 4 ) << " (avg.: " << ibf_config.true_avg_fp << ")" << std::endl;
    // 打印过滤器大小，以MB为单位，精确到2位小数
    std::cerr << std::fixed << std::setprecision( 2 ) << " - filter size: "
              << ( optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits ) / static_cast< double >( 8388608u )
              << "MB" << std::endl;
}

void print_stats_verbose( const StopClock& timeGanonBuild,
                          const StopClock& timeCountStoreHashes,
                          const StopClock& timeEstimateParams,
                          const StopClock& timeBuildIBF,
                          const StopClock& timeWriteIBF )
{
    /*
     * 打印构建过程中各个阶段的详细统计信息和时间
     */
    using ::operator<<;
    // 打印计算和存储哈希值的开始时间
    std::cerr << "Count/save hashes start: " << StopClock_datetime( timeCountStoreHashes.begin() ) << std::endl;
    // 打印计算和存储哈希值的结束时间
    std::cerr << "                    end: " << StopClock_datetime( timeCountStoreHashes.end() ) << std::endl;
    // 打印计算和存储哈希值的耗时（秒）
    std::cerr << "            elapsed (s): " << timeCountStoreHashes.elapsed() << std::endl;
    // 打印估计参数的开始时间
    std::cerr << "Estimate params   start: " << StopClock_datetime( timeEstimateParams.begin() ) << std::endl;
    // 打印估计参数的结束时间
    std::cerr << "                    end: " << StopClock_datetime( timeEstimateParams.end() ) << std::endl;
    // 打印估计参数的耗时（秒）
    std::cerr << "            elapsed (s): " << timeEstimateParams.elapsed() << std::endl;
    // 打印构建过滤器的开始时间
    std::cerr << "Building filter   start: " << StopClock_datetime( timeBuildIBF.begin() ) << std::endl;
	// 打印构建过滤器的结束时间
    std::cerr << "                    end: " << StopClock_datetime( timeBuildIBF.end() ) << std::endl;
	// 打印构建过滤器的耗时（秒）
    std::cerr << "            elapsed (s): " << timeBuildIBF.elapsed() << std::endl;
    // 打印保存过滤器的开始时间
    std::cerr << "Saving filer      start: " << StopClock_datetime( timeWriteIBF.begin() ) << std::endl;
    // 打印保存过滤器的结束时间
    std::cerr << "                    end: " << StopClock_datetime( timeWriteIBF.end() ) << std::endl;
    // 打印保存过滤器的耗时（秒）
    std::cerr << "            elapsed (s): " << timeWriteIBF.elapsed() << std::endl;
    // 打印整个构建过程的开始时间
    std::cerr << "ganon-build       start: " << StopClock_datetime( timeGanonBuild.begin() ) << std::endl;
    // 打印整个构建过程的结束时间
    std::cerr << "                    end: " << StopClock_datetime( timeGanonBuild.end() ) << std::endl;
    // 打印整个构建过程的耗时（秒）
    std::cerr << "            elapsed (s): " << timeGanonBuild.elapsed() << std::endl;


    std::cerr << std::endl;
}

} // namespace detail

bool run( Config config )
{

    // 验证用户提供的参数
    if ( !config.validate() )
        return false;

    // 如果启用了详细输出，则打印参数
    if ( config.verbose )
        std::cerr << config;

    // 开始计时器以记录总构建时间
    StopClock timeGanonBuild;
    StopClock timeCountStoreHashes;
    StopClock timeEstimateParams;
    StopClock timeBuildIBF;
    StopClock timeWriteIBF;
    timeGanonBuild.start();

    // 初始化统计信息（通用）和总计（每个线程一个）
    detail::Stats                stats;
    std::vector< detail::Total > totals( config.threads );

    // 创建 IBF 配置并设置固定参数
    IBFConfig ibf_config;
    ibf_config.kmer_size   = config.kmer_size;
    ibf_config.window_size = config.window_size;

    // 映射以存储每个目标的哈希数 {target: count}
    detail::THashesCount hashes_count;

    // 解析有效的输入文件到一个队列（SafeQueue），队列中的元素为以目标为键的输入文件映射
    SafeQueue< detail::InputFileMap > ifm_queue;
    ;
    // 调用 detail::parse_input_file 函数，解析输入文件，并初始化每个目标或序列ID的哈希计数
    for ( auto const& [target, files] :
          detail::parse_input_file( config.input_file, hashes_count, config.quiet, stats ) )
    {
        // 将目标和文件添加到 SafeQueue
        ifm_queue.push( detail::InputFileMap{ target, files } );
    }
    // 通知队列，推送操作已完成
    ifm_queue.notify_push_over();

    // 如果队列为空，表示没有有效的输入文件
    if ( ifm_queue.empty() )
    {
        std::cerr << "No valid input files" << std::endl;
        return false;
    }

    // 如果配置中指定了临时输出文件夹且该文件夹不存在，则创建它
    if ( config.tmp_output_folder != "" && !std::filesystem::exists( config.tmp_output_folder ) )
    {
        std::filesystem::create_directory( config.tmp_output_folder );
    }
    else
    {
        // 如果临时输出文件夹已存在，则删除之前创建的 .min 哈希文件
        detail::delete_hashes( hashes_count, config.tmp_output_folder );
    }

    // 开始计时，用于哈希计数和存储
    timeCountStoreHashes.start();
    // 初始化一个向量，用于存储并行任务的 future 对象
    std::vector< std::future< void > > tasks_count;

    // 创建并行任务来处理哈希计数
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        // 使用 std::async 启动并行任务，调用 detail::count_hashes 函数
        tasks_count.emplace_back( std::async( std::launch::async,
                                              detail::count_hashes,
                                              std::ref( ifm_queue ),
                                              std::ref( hashes_count ),
                                              std::ref( config ),
                                              std::ref( totals[taskn] ) ) );
    }
    // 等待所有线程完成
    for ( auto&& task : tasks_count )
    {
        // 获取每个任务的结果，确保任务完成
        task.get();
    }
    // 停止哈希计数和存储的计时
    timeCountStoreHashes.stop();
    // 汇总所有线程的总计数
    stats.add_totals( totals );

    // 根据提供的过滤器大小或最大假阳性率定义最优参数
    // 填充 ibf_config 结构体中的最优值
    timeEstimateParams.start();
    // 计算每个 bin 的最优最大哈希数（基于过滤器大小和fp之间的最小调和平均数）
    detail::optimal_hashes( config.max_fp,
                            config.filter_size,
                            ibf_config,
                            hashes_count,
                            config.hash_functions,
                            config.max_hash_functions,
                            config.mode );
    // 根据选择的 ibf_config 参数计算真实的最终假阳性率和平均假阳性率
    std::tie( ibf_config.true_max_fp, ibf_config.true_avg_fp ) = detail::true_false_positive(
        hashes_count, ibf_config.max_hashes_bin, ibf_config.bin_size_bits, ibf_config.hash_functions );
    // 停止参数估算的计时
    timeEstimateParams.stop();

    // 如果启用了详细模式，则打印IBF配置信息
    if ( config.verbose )
    {
        // 输出IBF配置的详细信息
        std::cerr << ibf_config;
        // 输出过滤器大小，以位和兆字节为单位
        std::cerr << "Filter size: " << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                  << " Bits";
        std::cerr << " ("
                  << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                         / static_cast< double >( 8388608u )
                  << " Megabytes)" << std::endl;
    }

    // 如果没有有效的bin数量，输出错误信息并返回false
    if ( ibf_config.n_bins == 0 )
    {
        std::cerr << "No valid sequences to build" << std::endl;
        return false;
    }


    // 将哈希值拆分成最佳大小的技术bin
    // bin_map_hash 映射 bin 编号到目标、起始哈希索引和结束哈希索引的元组
    const detail::TBinMapHash bin_map_hash = detail::create_bin_map_hash( ibf_config, hashes_count );

    // 创建IBF对象并开始计时
    timeBuildIBF.start();
    auto ibf = detail::TIBF{ seqan3::bin_count{ ibf_config.n_bins },
                             seqan3::bin_size{ ibf_config.bin_size_bits },
                             seqan3::hash_function_count{ ibf_config.hash_functions } };

    // 将最小化哈希值添加到IBF中，读取之前写入的.min文件
    // 并行化处理每批次（64个）以保证IBF的线程安全
    uint64_t batch_size = 64;   // 定义每批次处理的大小为64
    // 计算处理所有bin所需的最大批次数
    uint64_t max_batch = std::ceil( bin_map_hash.size() / static_cast< double >( batch_size ) );
    // 使用原子变量 bin_batches 跟踪线程中的批处理进度
    std::atomic< std::size_t >         bin_batches = 0;
    std::vector< std::future< void > > tasks_build;
    // 为每个线程创建并行任务来构建IBF
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        tasks_build.emplace_back( std::async( std::launch::async, [=, &ibf, &bin_batches, &bin_map_hash]() {
            detail::build( ibf, bin_batches, max_batch, batch_size, bin_map_hash, config.tmp_output_folder );
        } ) );
    }
    // 等待所有线程完成
    for ( auto&& task : tasks_build )
    {
        task.get();
    }
    timeBuildIBF.stop();

    // 删除临时的.min哈希文件
    detail::delete_hashes( hashes_count, config.tmp_output_folder );

    // 将IBF和其他数据结构写入文件
    timeWriteIBF.start();
    detail::save_filter( config, ibf, ibf_config, hashes_count, bin_map_hash );
    timeWriteIBF.stop(); // 停止写入操作的计时


    // 停止总构建时间的计时
    timeGanonBuild.stop();

    // 打印统计信息和时间
    if ( !config.quiet )
    {
        if ( config.verbose )
        {
            detail::print_stats_verbose(
                timeGanonBuild, timeCountStoreHashes, timeEstimateParams, timeBuildIBF, timeWriteIBF );
        }
        detail::print_stats( stats, ibf_config, timeGanonBuild );
    }

    return true;
}

} 
