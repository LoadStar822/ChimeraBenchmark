# Benchmark 数据与结果说明

本文件说明 ChimeraBenchmark 的数据、参考库、评估任务、指标和结果表。
它面向读完论文后希望进一步查看 benchmark 细节的读者。

## 目录

- [数据与结果导航](#数据与结果导航)
- [结果使用说明](#结果使用说明)
- [Reference Database Summary](#reference-database-summary)
- [Benchmark Dataset Summary](#benchmark-dataset-summary)
- [评估任务与口径](#评估任务与口径)
- [指标说明](#指标说明)
- [工具说明](#工具说明)
- [软件版本](#软件版本)

## 数据与结果导航

- `results/builds/README.md`：数据库构建成本。
- `results/classify/README.md`：逐读段分类准确性。
- `results/profile/README.md`：样本丰度画像准确性。
- `resources/reports/db_catalog.tsv`：Reference Database Summary 的 TSV 版本。
- `resources/reports/dataset_catalog.tsv`：Benchmark Dataset Summary 的 TSV 版本。

默认报告的分类层级为 `species` 和 `genus`。

## 结果使用说明

- 所有工具使用相同线程数；当前公开表格按 `threads=32` 解释。
- 除必要的输出格式选项外，分类和丰度推断使用工具默认参数。
- 评估统一使用 NCBI 分类体系（NCBI taxonomy）快照 `resources/taxonomy/ncbi_20260408/`。
- `ganon` 在结果表中显示为 `ganon2`，对应 ganon2 软件版本。
- Taxor 当前结果来自早期运行；引用主表结论前建议使用同一设置重新运行。

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
统计来自 `configs/datasets/*.yaml` 和输入序列文件扫描。

| Dataset Name | Total Size (GB) | Samples | Input Type | Reads / Contigs | Base Pairs (bp) | Mean Length (bp) | N50 (bp) | GC (%) | Q30 (%) | Truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAMI II Marine (long-read contigs, sample 0) | 0.9 | 1 | contig FASTA | 133,050 | 925,503,357 | 6,956 | 10,524 | 49.19 | — | per-read mapping + profile |
| CAMI II Marine (long-read contigs) | 10.8 | 10 | contig FASTA | 1,622,375 | 10,587,029,583 | 6,526 | 9,364 | 48.12 | — | per-read mapping + profile |
| CAMI II Marine (short-read contigs) | 10.8 | 10 | contig FASTA | 18,057,427 | 10,389,109,998 | 575 | 770 | 48.11 | — | per-read mapping + profile |
| ATCC MSA-1003 (Illumina) | 3.3 | 1 | paired FASTQ | 10,038,314 | 1,254,789,250 | 125 | 125 | 48.91 | 92.54 | profile only |
| ATCC MSA-1003 (PacBio HiFi) | 41.3 | 1 | single FASTQ | 2,419,037 | 20,543,987,923 | 8,493 | 9,198 | 53.75 | 96.70 | profile only |
| ZymoBIOMICS Even (GridION) | 28.5 | 1 | single FASTQ | 3,491,078 | 14,007,156,825 | 4,012 | 5,213 | 46.02 | 0.00 | profile only |
| ZymoBIOMICS Log (GridION) | 32.6 | 1 | single FASTQ | 3,667,007 | 16,032,264,247 | 4,372 | 5,290 | 40.94 | 0.00 | profile only |
| ZymoBIOMICS Even (PromethION) | 298.4 | 1 | single FASTQ | 36,527,376 | 146,291,709,850 | 4,005 | 5,288 | 46.38 | 0.00 | profile only |
| ZymoBIOMICS Log (PromethION) | 301.6 | 1 | single FASTQ | 35,118,078 | 148,028,914,788 | 4,215 | 5,219 | 41.41 | 0.00 | profile only |

## 评估任务与口径

本 benchmark 报告两个分类层级：`species` 和 `genus`。
逐读段分类（per-read classification）评估每条 read 或 contig 的分类结果。
丰度画像（abundance profiling）评估整个样本的物种组成估计。

在 `species` 或 `genus` 层级计算逐读段结果时，先将真值分类编号（taxid）和预测 taxid 都提升到同一目标层级，再判断是否一致。
这样可以公平处理原始真值为 strain、subspecies 或 isolate 的情况。

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

### Centrifuger

Centrifuger 的逐读段分类使用默认推断参数。
Profile 结果由 `centrifuger-quant` 从 classify 输出生成，并使用 CAMI profile 格式导出。

### Taxor

Taxor 结果来自早期运行，当前表格保留用于追踪已有实验。
在引用论文主表结论时，建议使用重新运行后的 Taxor 结果。

## 软件版本

- kraken2: 2.1.3
- bracken: 2.9
- centrifuger: 1.1.0-r291
- ganon2: 2.1.0
- sylph: 0.8.1
- taxor: 0.1.3（SeqAn 3.4.0-rc.1）
