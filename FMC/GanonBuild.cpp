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


// ����δѹ�������ݲ��ֵĽ���¡����������
typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > TIBF;

// �����ϣ����������ӳ������
typedef robin_hood::unordered_map< std::string, uint64_t > THashesCount;

// ����binӳ���ϣ������ӳ������
// �����ӳ���uint64_t���͵Ĺ�ϣֵ��һ������Ŀ���ַ���������uint64_t����������Ԫ��
typedef robin_hood::unordered_map< uint64_t, std::tuple< std::string, uint64_t, uint64_t > > TBinMapHash;

// ����һ���ṹ������ʾ�����ļ���ӳ��
struct InputFileMap
{
    std::string                target;  // Ŀ������
    std::vector< std::string > files;   // �ļ��б�
};



// ����һ���ṹ������ʾͳ�Ƶ�����
struct Total
{
    uint64_t files             = 0; // �ļ���
    uint64_t invalid_files     = 0; // ��Ч�ļ���
    uint64_t sequences         = 0; // ������
    uint64_t skipped_sequences = 0; // ������������
    uint64_t length_bp         = 0; // ���г��ȣ��Լ����Ϊ��λ��
};

// ����һ���ṹ������ʾͳ����Ϣ
struct Stats
{
    Total total;    // ����
    void  add_totals( std::vector< Total > const& totals )
    {
        // ������̵߳�������ӵ�ͳ����Ϣ��
        for ( auto const& t : totals )
        {
            total.files += t.files;     // �ۼ��ļ���
            total.invalid_files += t.invalid_files;     // �ۼ���Ч�ļ���
            total.sequences += t.sequences;     // �ۼ�������
            total.skipped_sequences += t.skipped_sequences;     // �ۼ�������������
            total.length_bp += t.length_bp;     // �ۼ����г���
        }
    }
};

inline std::string get_seqid( std::string header )
{
    /*
     * ���� FASTA ��ͷ�е����� ID��ֱ����һ���ո�֮ǰ���������ݣ�
     */
    return header.substr( 0, header.find( ' ' ) );
}

robin_hood::unordered_map< std::string, std::vector< std::string > > parse_input_file( const std::string& input_file,
                                                                                       THashesCount&      hashes_count,
                                                                                       bool               quiet,
                                                                                       Stats&             stats )
{
    /*
     * ���������ļ��ĺ��� -> һ�������ֶεı���ļ�: �ļ� [<tab> Ŀ�� <tab> ���� ID]
     * ����һ��ӳ�� {Ŀ��: {�ļ�: [���� ID]}}
     * ����ṩ�˵����е����� ID���� seqids ��һ������ ID
     * ������ļ���������ÿ���ļ�ֻ��һ��Ŀ�꣬seqids Ϊ��
     */
    robin_hood::unordered_map< std::string, std::vector< std::string > > input_map;
    std::string                                                          line;
    robin_hood::unordered_set< std::string >                             files;
    std::ifstream                                                        infile( input_file );
    // ���ж�ȡ�����ļ�
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        // ʹ���Ʊ���ָ�ÿ���ֶ�
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // file [<tab> Ŀ�� <tab> ���� ID]

        const std::string file = fields[0];
        files.insert( file );
        // ����ļ��Ƿ���ڻ��ļ��Ƿ�Ϊ��
        if ( !std::filesystem::exists( file ) || std::filesystem::file_size( file ) == 0 )
        {
            if ( !quiet )
                std::cerr << "WARNING: input file not found/empty: " << file << std::endl;
            stats.total.invalid_files++;
            continue;
        }

        if ( fields.size() == 1 )
        {
            // target ���ļ��������ļ�����û��·����
            const auto target = std::filesystem::path( file ).filename();
            input_map[target].push_back( file );
            hashes_count[target] = 0;
        }
        else if ( fields.size() == 2 )
        {
            // �ļ����ṩ�� target
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
     * �������еĹ�ϣֵ�洢��ָ���ļ����еĴ�����
     * ���δָ���ļ��У���洢�ڵ�ǰ�ļ���(".")��
     */
    // ��������ļ�·��
    std::filesystem::path outf{ tmp_output_folder };
    outf += target + ".min";    // �����ļ���Ϊ "target.min"
    // �Զ����ƺ�׷��ģʽ������ļ�
    std::ofstream outfile{ outf, std::ios::binary | std::ios::app };
    // ��ÿ����ϣֵд���ļ�
    for ( auto&& h : hashes )
    {
        outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
    }
    // �ر��ļ�
    outfile.close();
}


std::vector< uint64_t > load_hashes( std::string file )
{
    /*
     * �Ӵ��̼��ع�ϣ�ļ�������һ��������Щ��ϣֵ������
     */
    uint64_t                hash;
    std::vector< uint64_t > hashes;
    std::ifstream           infile{ file, std::ios::binary };
    // ��ȡ�ļ��е����й�ϣֵ������������
    while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
        hashes.push_back( hash );
    return hashes;
}

void delete_hashes( const THashesCount& hashes_count, const std::string tmp_output_folder )
{
    /*
     * �Ӵ���ɾ����ϣ�ļ�
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
     * �ú����ڶ���߳��в������У���һ�� InputFileMap �����ж�ȡ����
     * ����ÿ��Ŀ�����С������ϣֵ������洢������
     * �����ᱣ����Щ�������Ա��һ��ȷ����ѹ�������С
     */


    // ÿ���߳�һ����ͼ�����ڼ�����С������ϣ
    auto minimiser_view = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( config.kmer_size ) } );

    while ( true )
    {
        // �ȴ�ֱ���ж�ȡ���ݿ��ã��������ͽ������Ҷ���Ϊ��
        InputFileMap ifm = ifm_queue.pop();

        // ����ڵ�����Ϊ�գ��˳��߳�
        if ( ifm.target == "" )
            break;

        robin_hood::unordered_set< uint64_t > target_hashes;

        // ���Ŀ��������ļ�
        for ( auto& file : ifm.files )
        {
            if ( target_hashes.size() >= 2000000 )
                break;
            try
            {
                // ���ļ�
                seqan3::sequence_file_input< raptor::dna4_traits,
                                             seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                    fin{ file };


                // ���ļ���ΪĿ�� - �ӿ��ܰ���������е��ļ����������й�ϣֵ
                // Ȼ����м����ʹ洢

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

                    // �ﵽ���޺�����ֹͣ����ǰ�ļ�
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
     * ���ṹ��� IBF ���浽�ļ��С��ṹ�屻ת��Ϊ��׼���ʽ��
     * �������ݰ�����
     * �汾��Ϣ:   tuple(major, minor, patch)
     * ibf_config: ���ڹ����������� IBFConfig �ṹ��
     * hashes_count: vector<tuple(string, int)> [(Ŀ��, ����)]
     * bin_map:     vector<tuple(int, string)> [(bin���, Ŀ��)]
     * ibf:         seqan3 �еĽ�֯��¡������
     */
    // ����һ������������ļ���
    std::ofstream               os( config.output_file, std::ios::binary );
    // ����һ���������������
    cereal::BinaryOutputArchive archive( os );

    // ����ÿ�� bin ��ӳ�䲢�洢�ڱ�׼��ṹ�� {bin���: Ŀ��}
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    for ( auto& [binno, vals] : bin_map_hash )
    {
        bin_map.push_back( { binno, std::get< 0 >( vals ) } );
    }

    // �� hashes_count �洢�ڱ�׼��ṹ�� {Ŀ��: ��ϣ����}
    std::vector< std::tuple< std::string, uint64_t > > hashes_count_std;
    for ( auto& [target, count] : hashes_count )
    {
        hashes_count_std.push_back( { target, count } );
    }

    // ���л���������
    archive( std::make_tuple( defaults::version_major, defaults::version_minor, defaults::version_patch ) );
    archive( ibf_config );
    archive( hashes_count_std );
    archive( bin_map );
    archive( ibf );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes )
{
    /*
     * �������������ʺ�Ҫ�����Ԫ���������� bin ��С (��¡������)��
     */
    return std::ceil( ( n_hashes * std::log( max_fp ) ) / std::log( 1.0 / std::pow( 2, std::log( 2 ) ) ) );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes, uint8_t hash_functions )
{
    /*
     * �������������ʡ�Ҫ�����Ԫ�������͹�ϣ�������������� bin ��С (��¡������)��
     */
    return std::ceil( n_hashes
                      * ( -hash_functions / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) ) );
}

uint8_t hash_functions_from_ratio( uint64_t bin_size_bits, uint64_t n_hashes )
{
    /*
     * ���� bin ��С��Ԫ������������ѹ�ϣ����������
     */
    return static_cast< uint8_t >( std::log( 2 ) * ( bin_size_bits / static_cast< double >( n_hashes ) ) );
}

uint8_t get_optimal_hash_functions( uint64_t bin_size_bits,
                                    uint64_t n_hashes,
                                    uint8_t  hash_functions,
                                    uint8_t  max_hash_functions )
{
    /*
     * �ú�������Ƿ�Ӧ�ø��ݱ��ʼ����ϣ�������� (hash_functions = 0)��
     * ��ȷ����ϣ��������������Χ�� (1-5)��
     */

    uint8_t optimal_hash_functions = hash_functions;
    // �����ϣ��������Ϊ0������ݱ��ʼ�����ѹ�ϣ��������
    if ( optimal_hash_functions == 0 )
        optimal_hash_functions = hash_functions_from_ratio( bin_size_bits, n_hashes );
    // ���������Ĺ�ϣ���������������ֵ��Ϊ0����������Ϊ���ֵ
    if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
        optimal_hash_functions = max_hash_functions;

    return optimal_hash_functions;
}


uint64_t number_of_bins( THashesCount const& hashes_count, uint64_t n_hashes )
{
    /*
     * ����Ԫ���������� IBF �� bin ������������ֵ� bin����
     */
    uint64_t n_bins = 0;
    // ����ÿ��Ŀ�꼰���Ӧ�Ĺ�ϣ����
    for ( auto const& [target, count] : hashes_count )
    {
        // ����Ŀ������� bin ���������ۼӵ��� bin ������
        n_bins += std::ceil( count / static_cast< double >( n_hashes ) );
    }
    return n_bins;
}


double correction_rate( uint64_t max_split_bins, double max_fp, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * ���� bin �Ĵ�СӦ���ӵı��ʣ�����Ӧ�ɽ�Ŀ���ֳɶ�� bin �������Ķ��ش������⡣
     * ����Ŀ�������� bin �в�ֵ������
     */

    // ����Ŀ�������
    double const target_fpr        = 1.0 - std::exp( std::log( 1.0 - max_fp ) / max_split_bins );
    // �����µ� bin ��С
    size_t const new_bin_size      = bin_size( target_fpr, n_hashes, hash_functions );
    // ����ԭʼ�� bin ��С
    size_t const original_bin_size = bin_size( max_fp, n_hashes, hash_functions );
    // �����µ� bin ��С��ԭʼ bin ��С�ı���
    return static_cast< double >( new_bin_size ) / original_bin_size;
}


inline uint64_t optimal_bins( uint64_t n_bins )
{
    /*
     * ���ش��� IBF ������ bin ������64 �ı�����
     */
    return std::ceil( n_bins / 64.0 ) * 64;
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * ���ݲ������� bin ���������ʣ���¡��������
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
     * ���� IBF ���������ʣ�ƽ������󣩣�����targets��split�� bins
     */
    double highest_fp = 0;
    double average_fp = 0;
    // ����ÿ��Ŀ�������ʵ���ʣ����ǲ�ֳɶ�� bins����β��ԣ�
    for ( auto const& [target, count] : hashes_count )
    {
        // ʹ��ÿ�� bin ��ƽ����ϣ������������
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( max_hashes_bin ) );
        // ����ܻ���һ���ǳ�С�����������뵽��� bins��
        uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

        // ��ǰĿ������ʣ����ǲ�ֵ� bins
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
     * ���� `hashes_count` �е�����ϣ��
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
     * ����һ�����������ʻ��������С���������ܵ� bin ������������¡������������������Ӧ�Ĵ�С�����ǲ�ֵ� bin
     * ѡ������"���"֮����С�Ĺ�����/�������ʺ� bin �����Ĳ���������� ibf_config �ṹ
     */

    // ӵ�������С������Ŀ��
     uint64_t const max_hashes = get_max_hashes( hashes_count );

    // ����ƽ��ֵ����Сֵ
    uint64_t min_filter_size = 0;
    uint64_t min_bins        = 0;
    double   min_fp          = 1;

    // ��������ƽ��ֵ��ģ��ֵ
    struct SimParam
    {
        uint64_t n_hashes;
        uint64_t n_bins;
        uint64_t filter_size_bits;
        double   fp;
    };
    std::vector< SimParam > simulations;

    // ÿ 100 ��Ԫ�ؽ���һ��ģ��
    size_t iter = 100;
    // ��� max_hashes �Ƿ�С�ڵ�������
    if ( max_hashes < iter )
        iter = max_hashes;

    // (total + 1) ���ڴ���������
    for ( size_t n = max_hashes + 1; n > iter; n -= iter )
    {
        // Ҫ���뵽һ�� bin �е�Ԫ������
        uint64_t const n_hashes = n - 1;

        // ����Ŀ���Ԫ�ص�ʵ�� bin ���������� 64 �ı�����
        uint64_t const n_bins = number_of_bins( hashes_count, n_hashes );

        int64_t bin_size_bits          = 0;
        uint8_t optimal_hash_functions = 0;
        if ( filter_size )
        {
            // ����ṩ�� filter_size���򵥼���
            bin_size_bits = ( filter_size / static_cast< double >( optimal_bins( n_bins ) ) ) * 8388608u;
            optimal_hash_functions =
                get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
        }
        else
        {
            if ( hash_functions == 0 )
            {
                // ���ȶ����С��Ȼ�����ϣ��������
                bin_size_bits = bin_size( max_fp, n_hashes );
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
            }
            else
            {
                // �ṩ�˹�ϣ����������ʹ����������� bin ��С
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
                bin_size_bits = bin_size( max_fp, n_hashes, optimal_hash_functions );
            }
        }

        // Ŀ����Ϊ��� bin ��������
        uint64_t const max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        // ���� fp �� filter_size��ȡ�����ṩ�Ĳ���
        double   fp               = 0;
        uint64_t filter_size_bits = 0;
        if ( filter_size )
        {
            // ��Ե�ǰֵ����������ʣ����ǲ�ֵ� bin
            fp =
                1 - std::pow( 1.0 - false_positive( bin_size_bits, optimal_hash_functions, n_hashes ), max_split_bins );

            // save minimal value
            // ������Сֵ
            if ( fp < min_fp )
                min_fp = fp;
        }
        else
        {

            // ʵ�ʲ��뵽ÿ�� bin ��Ԫ������
            uint64_t const avg_n_hashes = std::ceil( max_hashes / static_cast< double >( max_split_bins ) );
            // ����ÿ����� bin ��ƽ�� n_hashes ������ʵ��������
            // ��������ã����߹������ʣ���� bin û����ȫ������
            double approx_fp = false_positive( bin_size_bits, optimal_hash_functions, avg_n_hashes );
            // �������ֵ�ϸߣ����ȼ��㣩�������
            if ( approx_fp > max_fp )
                approx_fp = max_fp;

            // ���ڵ���Ŀ�������������������
            double const crate = correction_rate( max_split_bins, approx_fp, optimal_hash_functions, n_hashes );
            // Ӧ���� bin ��С
            bin_size_bits = bin_size_bits * crate;
            // �������չ�������С
            filter_size_bits = bin_size_bits * optimal_bins( n_bins );

            // ֵ��С������С�� n_hashes ����ߵ� crate ����ֵ����ʱ�˳�ѭ��
            if ( filter_size_bits == 0 || std::isinf( crate ) )
                break;

            // ������Сֵ
            if ( filter_size_bits < min_filter_size || min_filter_size == 0 )
                min_filter_size = filter_size_bits;
        }

        // ����ģ��ֵ
        simulations.emplace_back( SimParam{ n_hashes, n_bins, filter_size_bits, fp } );

        /*        std::cout << "n_hashes: " << n_hashes << '\t';
                std::cout << "n_bins: " << n_bins << '\t';
                std::cout << "filter_size_bits: " << filter_size_bits << '\t';
                std::cout << "fp: " << fp << '\t';
                std::cout << "hash_functions: " << unsigned( optimal_hash_functions ) << '\n';*/

        // ������Сֵ
        if ( n_bins < min_bins || min_bins == 0 )
            min_bins = n_bins;
    }

    // ѡ�����š���ϣ��Ϊ n_bins �� filter_size �ĵ���ƽ����
    // ������������ܵ���Сֵ�Ĳ���
    // ���ѡ��������ģʽ��ƽ��ֵ�ᱻһ������ƫ�ƣ����ڱ�����˵��С���ã����Ϊ 0.5��
    // 0 ��ʾ���Զ�������ʹ������ֵ��������С��
    double mode_val = 1; // 1 ��Ĭ��ģʽ avg��n_bins �� filter_size �ĵ���ƽ����
    if ( mode == "smaller" || mode == "faster" )
        mode_val = 0.5;
    else if ( mode == "smallest" || mode == "fastest" )
        mode_val = 0;

    // ���� fp �� filter_size_bits��ȡ�����ṩ�Ĳ���
    double var_val = 1;
    // ���� bin ����
    double bins_val = 1;
    if ( mode == "smaller" || mode == "smallest" )
        var_val = mode_val;
    else if ( mode == "faster" || mode == "fastest" )
        bins_val = mode_val;

    // ��ʼ����Сƽ��ֵΪ0
    double min_avg = 0;
    // �������е�ģ�����
    for ( auto const& params : simulations )
    {
        // �������ʣ����ڱȽϵ�ǰ��������Сֵ֮��ı���
        double var_ratio = 0;
        // �����Ƿ��ṩfilter_size�������������
        if ( filter_size )
            // ����ṩ��filter_size������㵱ǰ������������С�������ʵı���
            var_ratio = params.fp / min_fp;
        else
            // ���û���ṩfilter_size������㵱ǰ��������С����С��������С�ı���
            var_ratio = params.filter_size_bits / static_cast< double >( min_filter_size );

        // ����bins���ʣ���ǰbins��������Сbins�����ı���
        double const bins_ratio = params.n_bins / static_cast< double >( min_bins );
        // ����ƽ��ֵ�����ݸ���ģʽ��ֵ������Ȩ��
        // 1 + mode_val^2 ȷ��Ȩ�ز�Ϊ��
        // (var_ratio * bins_ratio) ���������
        // (var_val * var_ratio) + (bins_val * bins_ratio) �Ǽ�Ȩ����
        double const avg        = ( 1 + std::pow( mode_val, 2 ) )
                           * ( ( var_ratio * bins_ratio ) / ( ( var_val * var_ratio ) + ( bins_val * bins_ratio ) ) );


        // �����ǰƽ��ֵС����Сƽ��ֵ����Сƽ��ֵΪ0���������Сƽ��ֵ��ibf����
        if ( avg < min_avg || min_avg == 0 )
        {
            min_avg = avg;  // ������Сƽ��ֵ
            // �����ṩ�Ĳ�������ibf����
            if ( filter_size )
            {
                // ����ṩ��filter_size������㲢����ibf��bin��С������������
                ibf_config.bin_size_bits =
                    ( filter_size / static_cast< double >( optimal_bins( params.n_bins ) ) ) * 8388608u;
                ibf_config.max_fp = params.fp;
            }
            else
            {
                // ���û���ṩfilter_size����ʹ�õ�ǰ�����Ĺ�������С�͸���������������
                ibf_config.bin_size_bits = params.filter_size_bits / optimal_bins( params.n_bins );
                ibf_config.max_fp        = max_fp;
            }

            // ����ibf���õ�����ϣ����bins�����͹�ϣ��������
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
     * ����IBF�����͵�ǰ�Ĺ�ϣֵ����bin���
     * ����ϣֵ��ƽ��ֵ���䵽�ָ��bin��
     * ���һ������bin��ź͹�ϣֵ������ӳ�䣬���ڲ��뵽IBF��
     */
    uint64_t    binno = 0;  // ��ʼ��bin���
    TBinMapHash bin_map_hash;   // ����binӳ���ϣ��
    for ( auto const& [target, count] : hashes_count )
    {
        // ����Ŀ����Ҫ��ƽ��bin����
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        // ����ÿ��bin�а����Ĺ�ϣ����
        uint64_t n_hashes_bin  = std::ceil( count / static_cast< double >( n_bins_target ) );

        // ���������Ĺ�ϣ���������������ֵ��������Ϊ���ֵ
        if ( n_hashes_bin > ibf_config.max_hashes_bin )
            n_hashes_bin = ibf_config.max_hashes_bin;
        // ����ÿ��Ŀ���bin
        for ( uint64_t i = 0; i < n_bins_target; ++i )
        {
            uint64_t hashes_idx_st = i * n_hashes_bin;  // ���㵱ǰbin����ʼ����
            uint64_t hashes_idx_en = hashes_idx_st + n_hashes_bin - 1;  // ���㵱ǰbin�Ľ�������
            // �����ʼ���������ܹ�ϣ����������ѭ��
            if ( hashes_idx_st >= count )
                break;
            // ����������������ܹ�ϣ����������Ϊ�������
            if ( hashes_idx_en >= count )
                hashes_idx_en = count - 1;
            // ����ǰbin����Ϣ����ӳ���ϣ����
            bin_map_hash[binno] = std::make_tuple( target, hashes_idx_st, hashes_idx_en );
            binno++;    // ����bin���
        }
    }
    // ���Դ�����bin������Ԥ����ͬ
    assert( bin_map_hash.size() == ibf_config.n_bins );
    return bin_map_hash;    // ����binӳ���ϣ��
}

void build( TIBF&                       ibf,
            std::atomic< std::size_t >& bin_batches,
            const uint64_t              max_batch,
            const uint64_t              batch_size,
            const TBinMapHash&          bin_map_hash,
            std::string                 tmp_output_folder )
{
    /*
     * ����С�����ļ�����ǰ���ɵ�ӳ�乹�� IBF
     */
    while ( true )
    {
        // ����ԭ�Ӽ����� bin_batches ���洢��ǰ���α��
        uint64_t batch = bin_batches++;
        if ( batch >= max_batch )
            break;

        // ���ò�������εı߽�
        uint64_t batch_start = batch * batch_size;
        uint64_t batch_end   = batch_start + batch_size - 1;
        if ( batch_end > bin_map_hash.size() - 1 )
            batch_end = bin_map_hash.size() - 1;

        // �洢�ļ��͹�ϣ��һ��ӳ�����Ա���ٷ���
        // ��Ϊ bin ���ο��ܻ����ظ��Ͷ���ļ�
        robin_hood::unordered_map< std::string, std::vector< uint64_t > > target_hashes;

        // ����������ϣ���뵽 ibf ��
        for ( uint64_t binno = batch_start; binno <= batch_end; binno++ )
        {
            auto [target, hashes_idx_st, hashes_idx_en] = bin_map_hash.at( binno );
            // ֻ��ȡһ���ļ����洢
            if ( !target_hashes.count( target ) )
            {
                auto file             = tmp_output_folder + target + ".min";
                target_hashes[target] = load_hashes( file );
            }
            // ����ϣ���뵽��Ӧ�� bin ��
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
     * ��ӡ����������ͳ����Ϣ
     */
    double elapsed = timeGanonBuild.elapsed();
    std::cerr << "ganon-build processed " << stats.total.sequences << " sequences / " << stats.total.files << " files ("
              << stats.total.length_bp / 1000000.0 << " Mbp) in " << elapsed << " seconds ("
              << ( stats.total.length_bp / 1000000.0 ) / ( elapsed / 60.0 ) << " Mbp/m)" << std::endl;

    // ���������Ч�ļ�����ӡ��Ч�ļ�����
    if ( stats.total.invalid_files > 0 )
        std::cerr << " - " << stats.total.invalid_files << " invalid files skipped" << std::endl;
    // ����������������У���ӡ��������������
    if ( stats.total.skipped_sequences > 0 )
        std::cerr << " - " << stats.total.skipped_sequences << " sequences skipped" << std::endl;

    // ��ӡ����ƽ���ļ������ʣ���ȷ��4λС��
    std::cerr << std::fixed << std::setprecision( 4 ) << " - max. false positive: " << ibf_config.true_max_fp;
    std::cerr << std::fixed << std::setprecision( 4 ) << " (avg.: " << ibf_config.true_avg_fp << ")" << std::endl;
    // ��ӡ��������С����MBΪ��λ����ȷ��2λС��
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
     * ��ӡ���������и����׶ε���ϸͳ����Ϣ��ʱ��
     */
    using ::operator<<;
    // ��ӡ����ʹ洢��ϣֵ�Ŀ�ʼʱ��
    std::cerr << "Count/save hashes start: " << StopClock_datetime( timeCountStoreHashes.begin() ) << std::endl;
    // ��ӡ����ʹ洢��ϣֵ�Ľ���ʱ��
    std::cerr << "                    end: " << StopClock_datetime( timeCountStoreHashes.end() ) << std::endl;
    // ��ӡ����ʹ洢��ϣֵ�ĺ�ʱ���룩
    std::cerr << "            elapsed (s): " << timeCountStoreHashes.elapsed() << std::endl;
    // ��ӡ���Ʋ����Ŀ�ʼʱ��
    std::cerr << "Estimate params   start: " << StopClock_datetime( timeEstimateParams.begin() ) << std::endl;
    // ��ӡ���Ʋ����Ľ���ʱ��
    std::cerr << "                    end: " << StopClock_datetime( timeEstimateParams.end() ) << std::endl;
    // ��ӡ���Ʋ����ĺ�ʱ���룩
    std::cerr << "            elapsed (s): " << timeEstimateParams.elapsed() << std::endl;
    // ��ӡ�����������Ŀ�ʼʱ��
    std::cerr << "Building filter   start: " << StopClock_datetime( timeBuildIBF.begin() ) << std::endl;
	// ��ӡ�����������Ľ���ʱ��
    std::cerr << "                    end: " << StopClock_datetime( timeBuildIBF.end() ) << std::endl;
	// ��ӡ�����������ĺ�ʱ���룩
    std::cerr << "            elapsed (s): " << timeBuildIBF.elapsed() << std::endl;
    // ��ӡ����������Ŀ�ʼʱ��
    std::cerr << "Saving filer      start: " << StopClock_datetime( timeWriteIBF.begin() ) << std::endl;
    // ��ӡ����������Ľ���ʱ��
    std::cerr << "                    end: " << StopClock_datetime( timeWriteIBF.end() ) << std::endl;
    // ��ӡ����������ĺ�ʱ���룩
    std::cerr << "            elapsed (s): " << timeWriteIBF.elapsed() << std::endl;
    // ��ӡ�����������̵Ŀ�ʼʱ��
    std::cerr << "ganon-build       start: " << StopClock_datetime( timeGanonBuild.begin() ) << std::endl;
    // ��ӡ�����������̵Ľ���ʱ��
    std::cerr << "                    end: " << StopClock_datetime( timeGanonBuild.end() ) << std::endl;
    // ��ӡ�����������̵ĺ�ʱ���룩
    std::cerr << "            elapsed (s): " << timeGanonBuild.elapsed() << std::endl;


    std::cerr << std::endl;
}

} // namespace detail

bool run( Config config )
{

    // ��֤�û��ṩ�Ĳ���
    if ( !config.validate() )
        return false;

    // �����������ϸ��������ӡ����
    if ( config.verbose )
        std::cerr << config;

    // ��ʼ��ʱ���Լ�¼�ܹ���ʱ��
    StopClock timeGanonBuild;
    StopClock timeCountStoreHashes;
    StopClock timeEstimateParams;
    StopClock timeBuildIBF;
    StopClock timeWriteIBF;
    timeGanonBuild.start();

    // ��ʼ��ͳ����Ϣ��ͨ�ã����ܼƣ�ÿ���߳�һ����
    detail::Stats                stats;
    std::vector< detail::Total > totals( config.threads );

    // ���� IBF ���ò����ù̶�����
    IBFConfig ibf_config;
    ibf_config.kmer_size   = config.kmer_size;
    ibf_config.window_size = config.window_size;

    // ӳ���Դ洢ÿ��Ŀ��Ĺ�ϣ�� {target: count}
    detail::THashesCount hashes_count;

    // ������Ч�������ļ���һ�����У�SafeQueue���������е�Ԫ��Ϊ��Ŀ��Ϊ���������ļ�ӳ��
    SafeQueue< detail::InputFileMap > ifm_queue;
    ;
    // ���� detail::parse_input_file ���������������ļ�������ʼ��ÿ��Ŀ�������ID�Ĺ�ϣ����
    for ( auto const& [target, files] :
          detail::parse_input_file( config.input_file, hashes_count, config.quiet, stats ) )
    {
        // ��Ŀ����ļ���ӵ� SafeQueue
        ifm_queue.push( detail::InputFileMap{ target, files } );
    }
    // ֪ͨ���У����Ͳ��������
    ifm_queue.notify_push_over();

    // �������Ϊ�գ���ʾû����Ч�������ļ�
    if ( ifm_queue.empty() )
    {
        std::cerr << "No valid input files" << std::endl;
        return false;
    }

    // ���������ָ������ʱ����ļ����Ҹ��ļ��в����ڣ��򴴽���
    if ( config.tmp_output_folder != "" && !std::filesystem::exists( config.tmp_output_folder ) )
    {
        std::filesystem::create_directory( config.tmp_output_folder );
    }
    else
    {
        // �����ʱ����ļ����Ѵ��ڣ���ɾ��֮ǰ������ .min ��ϣ�ļ�
        detail::delete_hashes( hashes_count, config.tmp_output_folder );
    }

    // ��ʼ��ʱ�����ڹ�ϣ�����ʹ洢
    timeCountStoreHashes.start();
    // ��ʼ��һ�����������ڴ洢��������� future ����
    std::vector< std::future< void > > tasks_count;

    // �������������������ϣ����
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        // ʹ�� std::async �����������񣬵��� detail::count_hashes ����
        tasks_count.emplace_back( std::async( std::launch::async,
                                              detail::count_hashes,
                                              std::ref( ifm_queue ),
                                              std::ref( hashes_count ),
                                              std::ref( config ),
                                              std::ref( totals[taskn] ) ) );
    }
    // �ȴ������߳����
    for ( auto&& task : tasks_count )
    {
        // ��ȡÿ������Ľ����ȷ���������
        task.get();
    }
    // ֹͣ��ϣ�����ʹ洢�ļ�ʱ
    timeCountStoreHashes.stop();
    // ���������̵߳��ܼ���
    stats.add_totals( totals );

    // �����ṩ�Ĺ�������С�����������ʶ������Ų���
    // ��� ibf_config �ṹ���е�����ֵ
    timeEstimateParams.start();
    // ����ÿ�� bin ����������ϣ�������ڹ�������С��fp֮�����С����ƽ������
    detail::optimal_hashes( config.max_fp,
                            config.filter_size,
                            ibf_config,
                            hashes_count,
                            config.hash_functions,
                            config.max_hash_functions,
                            config.mode );
    // ����ѡ��� ibf_config ����������ʵ�����ռ������ʺ�ƽ����������
    std::tie( ibf_config.true_max_fp, ibf_config.true_avg_fp ) = detail::true_false_positive(
        hashes_count, ibf_config.max_hashes_bin, ibf_config.bin_size_bits, ibf_config.hash_functions );
    // ֹͣ��������ļ�ʱ
    timeEstimateParams.stop();

    // �����������ϸģʽ�����ӡIBF������Ϣ
    if ( config.verbose )
    {
        // ���IBF���õ���ϸ��Ϣ
        std::cerr << ibf_config;
        // �����������С����λ�����ֽ�Ϊ��λ
        std::cerr << "Filter size: " << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                  << " Bits";
        std::cerr << " ("
                  << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                         / static_cast< double >( 8388608u )
                  << " Megabytes)" << std::endl;
    }

    // ���û����Ч��bin���������������Ϣ������false
    if ( ibf_config.n_bins == 0 )
    {
        std::cerr << "No valid sequences to build" << std::endl;
        return false;
    }


    // ����ϣֵ��ֳ���Ѵ�С�ļ���bin
    // bin_map_hash ӳ�� bin ��ŵ�Ŀ�ꡢ��ʼ��ϣ�����ͽ�����ϣ������Ԫ��
    const detail::TBinMapHash bin_map_hash = detail::create_bin_map_hash( ibf_config, hashes_count );

    // ����IBF���󲢿�ʼ��ʱ
    timeBuildIBF.start();
    auto ibf = detail::TIBF{ seqan3::bin_count{ ibf_config.n_bins },
                             seqan3::bin_size{ ibf_config.bin_size_bits },
                             seqan3::hash_function_count{ ibf_config.hash_functions } };

    // ����С����ϣֵ��ӵ�IBF�У���ȡ֮ǰд���.min�ļ�
    // ���л�����ÿ���Σ�64�����Ա�֤IBF���̰߳�ȫ
    uint64_t batch_size = 64;   // ����ÿ���δ���Ĵ�СΪ64
    // ���㴦������bin��������������
    uint64_t max_batch = std::ceil( bin_map_hash.size() / static_cast< double >( batch_size ) );
    // ʹ��ԭ�ӱ��� bin_batches �����߳��е����������
    std::atomic< std::size_t >         bin_batches = 0;
    std::vector< std::future< void > > tasks_build;
    // Ϊÿ���̴߳�����������������IBF
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        tasks_build.emplace_back( std::async( std::launch::async, [=, &ibf, &bin_batches, &bin_map_hash]() {
            detail::build( ibf, bin_batches, max_batch, batch_size, bin_map_hash, config.tmp_output_folder );
        } ) );
    }
    // �ȴ������߳����
    for ( auto&& task : tasks_build )
    {
        task.get();
    }
    timeBuildIBF.stop();

    // ɾ����ʱ��.min��ϣ�ļ�
    detail::delete_hashes( hashes_count, config.tmp_output_folder );

    // ��IBF���������ݽṹд���ļ�
    timeWriteIBF.start();
    detail::save_filter( config, ibf, ibf_config, hashes_count, bin_map_hash );
    timeWriteIBF.stop(); // ֹͣд������ļ�ʱ


    // ֹͣ�ܹ���ʱ��ļ�ʱ
    timeGanonBuild.stop();

    // ��ӡͳ����Ϣ��ʱ��
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
