# Benchmark 数据与结果说明

本文件说明 ChimeraBenchmark 的数据、参考库、评估任务、指标和结果表。
它面向希望复核数据来源、评估口径和结果的读者。

## 目录

- [数据与结果导航](#数据与结果导航)
- [结果使用说明](#结果使用说明)
- [Reference Database Summary](#reference-database-summary)
- [Benchmark Dataset Summary](#benchmark-dataset-summary)
- [评估任务与口径](#评估任务与口径)
- [数据集特定口径](#数据集特定口径)
- [核心结果汇总](#核心结果汇总)
- [真实队列 clade-level 实验](#真实队列-clade-level-实验)
- [指标说明](#指标说明)
- [工具说明](#工具说明)
- [软件版本](#软件版本)

## 数据与结果导航

- `results/builds/README.md`：数据库构建成本。
- `results/classify/README.md`：逐读段分类准确性。
- `results/profile/README.md`：样本丰度画像准确性。
- `results/real/README.md`：真实队列 clade-level 实验结果。
- `results/paper_run_manifest.tsv`：汇总结果所对应的运行、工具版本和完成状态。
- `results/*/summary.tsv`：build、classify 和 profile 的机器可读汇总表。
- `resources/reports/db_catalog.tsv`：Reference Database Summary 的 TSV 版本。
- `resources/reports/dataset_catalog.tsv`：Benchmark Dataset Summary 的 TSV 版本。

默认报告的分类层级为 `species` 和 `genus`。

## 结果使用说明

- 默认 species/genus benchmark 中所有工具使用相同线程数，按 `threads=32` 解释；Fna C2 真实队列实验统一使用 128 threads，并在其 README 中单独记录。
- 除必要的输出格式选项外，分类和丰度推断使用工具默认参数。
- 评估统一使用 NCBI 分类体系（NCBI taxonomy）快照 `resources/taxonomy/ncbi_20260408/`。
- `ganon` 在结果表中显示为 `ganon2`，对应 ganon2 软件版本。
- Taxor 条目来自较早的运行批次；对应版本和设置见各数据集结果。

## Reference Database Summary

本表概括 benchmark 使用的参考库规模。
它描述参考库本体，而不是某个工具构建出的数据库目录。
统计来源为 reference target TSV 和 NCBI assembly summary 元数据。

| Dataset Name | Total Size (GB) | Total Sequences | Base Pairs (bp) | Assemblies | Species Count | Source |
| --- | --- | --- | --- | --- | --- | --- |
| CAMI RefSeq | 173.0 | 16,738,778 | 588,675,466,098 | 141,677 | 21,075 | NCBI RefSeq genomic, CAMI snapshot |
| RefSeq Complete | 55.3 | 119,431 | 189,355,217,616 | 58,666 | 25,126 | NCBI RefSeq complete genomes |

## Benchmark Dataset Summary

本表概括 benchmark 使用的实际输入文件。
统计来自 `configs/datasets/*.yaml`、输入序列文件扫描或数据集随附元数据。

| Dataset Name | Total Size (GB) | Samples | Input Type | Reads / Contigs | Base Pairs (bp) | Mean Length (bp) | N50 (bp) | GC (%) | Q30 (%) | Truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAMI II Marine (long-read contigs, sample 0) | 0.9 | 1 | contig FASTA | 133,050 | 925,503,357 | 6,956 | 10,524 | 49.19 | — | per-read mapping + profile |
| CAMI II Marine (long-read contigs) | 10.8 | 10 | contig FASTA | 1,622,375 | 10,587,029,583 | 6,526 | 9,364 | 48.12 | — | per-read mapping + profile |
| CAMI II Marine (short-read contigs) | 10.8 | 10 | contig FASTA | 18,057,427 | 10,389,109,998 | 575 | 770 | 48.11 | — | per-read mapping + profile |
| CAMI II Strain Madness (long reads, sample 0) | 4.0 | 1 | single FASTQ | 698,731 | 2,000,316,047 | 2,863 | 3,117 | 46.94 | 0.00 | per-read mapping; profile from read abundance |
| CAMI II Strain Madness (short reads, sample 0) | 4.2 | 1 | paired FASTQ | 13,310,046 | 1,996,506,900 | 150 | 150 | — | — | per-read mapping; profile from read abundance |
| CAMI II Strain Madness (long reads) | 401.1 | 100 | single FASTQ | 69,647,201 | — | — | — | — | — | per-read mapping; profile from read abundance |
| CAMI II Strain Madness (short reads) | 424.2 | 100 | paired FASTQ | 1,330,822,164 | 199,623,324,600 | 150 | 150 | — | — | per-read mapping; profile from read abundance |
| ATCC MSA-1003 (Illumina) | 3.3 | 1 | paired FASTQ | 10,038,314 | 1,254,789,250 | 125 | 125 | 48.91 | 92.54 | profile only |
| ATCC MSA-1003 (PacBio HiFi) | 41.3 | 1 | single FASTQ | 2,419,037 | 20,543,987,923 | 8,493 | 9,198 | 53.75 | 96.70 | profile only |
| ZymoBIOMICS Even (GridION) | 28.5 | 1 | single FASTQ | 3,491,078 | 14,007,156,825 | 4,012 | 5,213 | 46.02 | 0.00 | profile only |
| ZymoBIOMICS Log (GridION) | 32.6 | 1 | single FASTQ | 3,667,007 | 16,032,264,247 | 4,372 | 5,290 | 40.94 | 0.00 | profile only |
| ZymoBIOMICS Even (PromethION) | 298.4 | 1 | single FASTQ | 36,527,376 | 146,291,709,850 | 4,005 | 5,288 | 46.38 | 0.00 | profile only |
| ZymoBIOMICS Log (PromethION) | 301.6 | 1 | single FASTQ | 35,118,078 | 148,028,914,788 | 4,215 | 5,219 | 41.41 | 0.00 | profile only |
| PRJNA637878 Supported Fecal Genomes | 110.5 | 19 | paired FASTQ | 1,355,677,392 | 169,812,107,628 | 125 | — | — | — | local per-read + supported-strain profile |
| PRJNA637878 Supported Fecal Genomes (single-read) | 110.5 | 19 | single FASTQ | 1,321,826,777 | 169,812,107,628 | 128 | — | — | — | mate-level per-read + supported-strain profile |

## 评估任务与口径

本 benchmark 报告两个分类层级：`species` 和 `genus`。
逐读段分类（per-read classification）评估每条 read 或 contig 的分类结果。
丰度画像（abundance profiling）评估整个样本的物种组成估计。

在 `species` 或 `genus` 层级计算逐读段结果时，先将真值分类编号（taxid）和预测 taxid 都提升到同一目标层级，再判断是否一致。
这样可以公平处理原始真值为 strain、subspecies 或 isolate 的情况。

## 数据集特定口径

### PRJNA637878

PRJNA637878 的部分测序运行（SRA run）包含未配对读段（singleton reads），表现为 `spots_with_mates` 小于 `spots`。
`prjna637878-supported19-single-read` 将 R1 和 R2 作为独立单读段评估（single-read evaluation），并在生成 FASTQ 和真值表时为读段标识符（read identifier）增加 `__mate1` 或 `__mate2` 配对端标记（mate tag）。
该设置不改变原始实验的 paired-end 来源；它只规定评估时所有工具使用同一批单读段输入，且不使用 mate-pair 信息。

## 核心结果汇总

跨数据集汇总由输入和评估口径一致、状态为 `complete` 的运行生成；对应运行信息记录在 `results/paper_run_manifest.tsv`。
逐读段分类汇总覆盖 CAMI II Strain Madness long/short、CAMI II Marine long/short 和 `prjna637878-supported19-single-read`，比较 Chimera、Centrifuger 与 Kraken2。
丰度画像结果独立汇总工具原生输出的 sample-level abundance estimates。

机器可读结果位于 `results/builds/summary.tsv`、`results/classify/summary.tsv`、`results/classify/sample_metrics.tsv` 和 `results/profile/summary.tsv`。

## 真实队列 clade-level 实验

Fna C2 真实队列实验位于 `results/real/`。
该实验在默认 `species/genus` benchmark 之外，评估原论文定义的近邻 clade 信号及其 read-level evidence。

`results/real/README.md` 汇总三队列 Fna C2 分析的实验设计、指标、结果和配套机器表。

Detection floor 先根据原论文报告的 C2 abundance 和当前输入深度换算 expected C2 reads，再将达到指定下限的 C2-positive 样本纳入阳性组。
该口径用于区分“总体是否有 C2 信号”和“当 C2 信号达到可检测强度时，工具是否能把样本拉开”。

## 指标说明

### Per-read 指标

逐读段指标衡量单条序列是否被分到正确的目标分类层级。

- `Precision`：被工具分到某个目标层级的序列中，有多少是正确的。
- `Recall`：真值中可评估的序列里，有多少被正确找回。
- `F1`：precision 和 recall 的调和平均。
- `Truth Mapped Rate`：真值中有多少序列能映射到目标层级。
- `Pred Mapped Rate`：预测中有多少序列能映射到目标层级。

当 mapped rate 明显低于 `1` 时，F1 代表可映射子集上的表现，不应单独解释为全数据集表现。

### Profile 指标

Profile 表使用 OPAL 核心指标（OPAL core metrics）：
- `Completeness`：真值中出现的 taxa 有多少被预测到。
- `Purity`：预测到的 taxa 中有多少确实存在于真值。
- `L1 Norm`：真值和预测丰度分布之间的 L1 距离，范围为 `0..2`。
- `Weighted UniFrac`：沿分类体系树（taxonomy tree）比较两边的累积丰度差异。

Completeness 和 purity 越高越好；L1 Norm 和 Weighted UniFrac 越低越好。

Profile 结果只使用工具原生输出的丰度文件。
如果某工具没有原生 profile 输出，它不会出现在 profile 结果表中。

### 字段说明

- `Elapsed (s)`：运行耗时，单位为秒。
- `Max RSS (GB)`：最大常驻内存，单位为 GB。
- `Total Size (GB)`：输入文件或参考库文件的十进制 GB 大小。
- `Q30 (%)`：FASTQ 碱基质量分数不低于 30 的比例。

## 工具说明

### Bracken

Bracken 基于 Kraken2 输出进行丰度重估，不提供逐读段分类结果。
因此 Bracken 只出现在 profile 结果表中。
运行时间和内存按完整 `kraken2 -> bracken` 流程统计。

Bracken 的数据库需要额外生成 `database100mers.kmer_distrib`。
该步骤使用 Bracken 默认读长设置 `read_len=100`。

### ganon2

ganon2 在 CAMI II Strain Madness short reads 上使用固定批次（fixed-size batch）运行，当前为每批 2 个样本。
原因是默认期望最大化（expectation maximization, EM）重分配在全量拼接输入上超过可用内存。
该设置保留 ganon2 默认的多匹配处理方式，但 EM 的估计范围为每个 batch 内部。
结果表中这一路径仍对应 CAMI II Strain Madness short-read dataset；`Samples` 表示完成的 batch 数。

ganon2 不报告 PRJNA637878 single-read lane 结果。
该数据集在默认 EM 重分配阶段峰值内存过高，在当前服务器上未能完成运行。
失败样本的 Max RSS 约为 `495..596 GB`，进程在 reassigning reads 阶段被系统终止。

### Centrifuger

Centrifuger 的逐读段分类使用默认推断参数。
Profile 结果由 `centrifuger-quant` 从 classify 输出生成，并使用 CAMI profile 格式导出。

### Taxor

Taxor 条目来自早期运行，与当前跨工具汇总的运行批次不同。

## 软件版本

- kraken2: 2.1.3
- bracken: 2.9
- centrifuger: 1.1.0-r291
- ganon2: 2.1.0
- sylph: 0.8.1
- taxor: 0.1.3（SeqAn 3.4.0-rc.1）
