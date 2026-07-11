# Real-data Results

本目录保存默认 `species/genus` benchmark 之外的真实队列 clade-level 实验。
当前实验围绕已发表 CRC 队列中的 paper-defined Fna C2 信号，评估工具能否在近邻参考库下恢复该 clade 信号，并把样本层面的 abundance signal 转化为可审计的 read-level evidence。
结果分别报告 read classification、逐 read 审计和 sample-level signal，以对应不同工具的原生输出层级。

## Dataset: fna-c2-crc3-head3m

### Experiment Design

Fna C2 CRC stool metagenome 实验合并 Yachida、Yu 和 Wirbel 三个队列，共 760 个 R1 samples。
每个样本使用 R1 的 head 3M reads；源文件不足 3M reads 的样本使用可用 R1 reads。

样本标签来自原论文报告的 Fna C1/C2 abundance signal，而不是由 Chimera 或其他工具重新定义：

| Role | Samples | Meaning | Included in C2-positive vs zero-Fna AUC |
| --- | ---: | --- | --- |
| paper C2-positive | 137 | 原论文报告存在 Fna C2 signal | yes |
| paper zero-Fna | 619 | 原论文报告 Fna signal 为零 | yes |
| paper Fna-positive / C2-zero | 4 | 原论文报告有 Fna signal 但 C2 为零 | no |

参考库是冻结的 Fna near-neighbor panel，包含 paper-defined Fna C1、Fna C2、non-C1/C2 F. nucleatum comparators，以及其他 Fusobacterium near-neighbor / decoy references。

### Tool Settings

三个队列固定使用同一批输入 reads 和同一套 reference sequences。除线程数与输出路径外，分类或 profiling 阶段使用工具默认参数。

| Tool | Version | Run setting |
| --- | --- | --- |
| Chimera | 1.6.3, build 2026-06-27 | default `classify`, 128 threads; `--profile-cami` 只请求附加输出，不改变 read classification |
| Centrifuger | source `v1.0.12-10-gdf5bd32` | default single-read classification, 128 threads |
| Kraken2 | 2.1.3 | default classification, 128 threads; DB build 使用 load factor 0.1，因为 default 和 0.2 build 均触发 compact hash table capacity failure |
| sylph | 0.8.1 | default genome sketch parameter `c=200` and default `profile`, 128 threads |
| minimap2 audit | 2.24-r1122 | orthogonal short-read alignment audit |

### Metric Definitions

`expected C2 reads` 由原论文 Fna C2 abundance 按本实验输入深度换算：

```text
expected C2 reads = paper Fna C2 percent * input reads / 100
```

`detection floor` 指只把 `expected C2 reads` 达到指定下限的 paper C2-positive 样本作为阳性；zero-Fna 样本保持不变。
例如 `floor >=1000` 表示只评估在当前抽样深度下，原论文信号预计至少贡献 1000 条 C2 reads 的 positive samples。

`sample-level C2 signal` 的定义按工具输出类型区分：

| Tool type | Signal definition |
| --- | --- |
| Chimera / Centrifuger / Kraken2 LF01 | reads assigned to paper-defined Fna C2 per million input reads |
| sylph | sequence abundance percentage reported for paper-defined Fna C2 references |

`AUC` 衡量工具信号能否把 paper C2-positive 样本排在 zero-Fna 样本前面。
数值越高表示正负样本排序分离越清晰；`0.5` 接近随机排序，`1.0` 表示完整排序分离。

`Pearson vs paper C2` 是工具 C2 signal 与原论文 Fna C2 abundance signal 的线性相关性。

`read-level audit` 用正交 read alignment 审计工具叫出的 C2 candidate reads。
候选 reads 被重新比对到冻结的 paper-clade reference panel，并按最高分比对结果分为：

| Audit metric | Definition |
| --- | --- |
| raw candidate reads/M | 工具原始叫成 Fna C2 的 reads per million |
| audited C2-supported reads/M | 最高分比对包含 Fna C2 reference 的 candidate reads per million，包括 strict 和 tied support |
| audited strict-C2-best reads/M | 最高分比对只落在 Fna C2 reference 上的 candidate reads per million |

下列表格汇总三队列合并分析的主结果。

### Read-level Classification and Audit

Read-level audit 把工具叫出的 C2 candidate reads 投影到 alignment-supported C2 evidence。
同一个 AUC 口径被用于比较 positive samples 与 zero-Fna samples，只是信号从 raw calls 换成审计后的 reads/M。

#### Raw Candidate Reads per Million

| Tool | all positives | floor >=500 | floor >=1000 | floor >=2000 | floor >=5000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.613 | 0.807 | 0.881 | 0.945 | 0.938 |
| Centrifuger | 0.631 | 0.737 | 0.799 | 0.846 | 0.905 |
| Kraken2 LF01 | 0.532 | 0.604 | 0.651 | 0.704 | 0.762 |

#### Audited C2-supported Reads per Million

| Tool | all positives | floor >=500 | floor >=1000 | floor >=2000 | floor >=5000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.830 | 0.966 | 0.981 | 0.992 | 0.995 |
| Centrifuger | 0.661 | 0.804 | 0.873 | 0.921 | 0.959 |
| Kraken2 LF01 | 0.530 | 0.601 | 0.647 | 0.700 | 0.757 |

#### Audited Strict-C2-best Reads per Million

| Tool | all positives | floor >=500 | floor >=1000 | floor >=2000 | floor >=5000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.881 | 0.988 | 0.995 | 0.998 | 0.999 |
| Centrifuger | 0.669 | 0.816 | 0.884 | 0.927 | 0.955 |
| Kraken2 LF01 | 0.527 | 0.597 | 0.642 | 0.695 | 0.751 |

10,000 次按 cohort 和 sample role 分层的 paired bootstrap 给出以下 strict-C2-best endpoint 置信区间。`Chimera minus tool` 使用每次相同的重采样样本计算，因此直接检验同一批样本上的 AUC 差值。

| Positive set | Tool | AUC (95% CI) | Chimera minus tool (95% CI) |
| --- | --- | ---: | ---: |
| all positives | Chimera | **0.881 (0.848 to 0.911)** | - |
| all positives | Centrifuger | 0.669 (0.618 to 0.718) | **0.212 (0.169 to 0.256)** |
| all positives | Kraken2 LF01 | 0.527 (0.478 to 0.576) | **0.353 (0.301 to 0.404)** |
| floor >=1000 | Chimera | **0.995 (0.990 to 0.998)** | - |
| floor >=1000 | Centrifuger | 0.884 (0.821 to 0.942) | **0.110 (0.055 to 0.170)** |
| floor >=1000 | Kraken2 LF01 | 0.642 (0.557 to 0.726) | **0.353 (0.270 to 0.437)** |

四个 paired differences 的 two-sided bootstrap P 均小于 `2e-4`。完整三个 read-level endpoints 和全部 detection floors 的置信区间位于对应 TSV。

### C2 Specificity

为直接检验 strict evidence 反映的是 C2 还是一般 Fna signal，在每个队列中同时用原论文 C2 和 C1 abundance 解释 `strict-C2-best reads / input reads`，并调整 CRC condition、sex、age、BMI 和 source sequencing depth。变量在队列内变换并标准化，随后合并 partial standardized beta。

| Published abundance predictor | Three-cohort partial beta | 95% CI | P value |
| --- | ---: | ---: | ---: |
| Fna C2 | **0.875** | **0.759 to 0.991** | **3.81e-49** |
| Fna C1 | 0.020 | -0.012 to 0.052 | 0.221 |

C2 partial beta 在 Yu、Wirbel 和 Yachida 分别为 `0.766`、`0.976` 和 `0.847`，三个队列均显著；同一模型中的 C1 beta 均不显著。这表明 Chimera 的 strict read evidence 跟随 paper-defined C2 abundance，而不是只跟随非特异的 Fna/C1 signal。

### CRC Association from Audited Read-level Evidence

为检验工具信号是否恢复原论文的 CRC 生物学关联，本分析直接比较全部 760 个样本中的 CRC 与 control，而不使用 paper C2-positive / zero-Fna 标签定义结局。
每个队列对 arcsine-square-root transformed endpoint 拟合与原论文一致的协变量模型：CRC condition、sex、age、BMI 和 source sequencing depth；随后用 Paule-Mandel random-effects model 合并三个队列的 adjusted standardized effect。

Read-level endpoint 为 `audited strict-C2-best reads/M`。
sylph 输出样本层面的丰度估计，不提供逐 read assignment，因此不适用于该 read-level audit。

| Tool / evidence | Three-cohort effect | 95% CI | P value | I2 |
| --- | ---: | ---: | ---: | ---: |
| Published Fna C2 abundance | 0.384 | 0.240 to 0.528 | 1.77e-7 | 0.0% |
| Chimera strict-C2-best reads/M | **0.336** | **0.193 to 0.480** | **4.56e-6** | **0.0%** |
| Centrifuger strict-C2-best reads/M | 0.064 | -0.179 to 0.307 | 0.604 | 45.1% |
| Kraken2 LF01 strict-C2-best reads/M | -0.275 | -0.780 to 0.230 | 0.286 | 86.0% |

Chimera 的 adjusted effect 在三个队列中方向一致，并逐队列达到显著：

| Cohort | Published Fna C2 effect | Chimera strict-read effect | Chimera P value |
| --- | ---: | ---: | ---: |
| Yu | 0.437 | **0.532** | 0.00402 |
| Wirbel | 0.438 | **0.406** | 0.0261 |
| Yachida | 0.358 | **0.273** | 0.00230 |

Secondary endpoint `strict-C2-best rate` 给出一致结果：Chimera three-cohort effect 为 `0.466`（95% CI `0.251 to 0.681`，P = `2.18e-5`）。
Chimera 因而把该生物学方向转化为逐 read 可审计的证据，并在三个队列中重建接近 published Fna C2 abundance 的 adjusted effect。

### Sample-level Signals

下列结果按各工具的原生输出统计 sample-level C2 signal。

| Tool | Signal unit | Samples | C2-positive mean | zero-Fna mean | all-positive AUC | Pearson vs paper C2 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Centrifuger | reads/M | 760 | 322.944 | 241.509 | 0.633 | 0.526 |
| Chimera | reads/M | 760 | 1360.832 | 978.246 | 0.613 | 0.853 |
| Kraken2 LF01 | reads/M | 760 | 79.601 | 70.806 | 0.532 | 0.260 |
| sylph | sequence abundance % | 760 | 15.299 | 0.000 | 0.613 | 0.587 |

| Tool | all positives (n=137) | floor >=500 (n=50) | floor >=1000 (n=38) | floor >=2000 (n=26) | floor >=5000 (n=12) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.613 | 0.807 | 0.881 | 0.945 | 0.938 |
| Centrifuger | 0.633 | 0.738 | 0.800 | 0.847 | 0.906 |
| sylph | 0.613 | 0.810 | 0.882 | 0.962 | 0.958 |
| Kraken2 LF01 | 0.532 | 0.604 | 0.651 | 0.704 | 0.762 |

Native sample-level CRC association 中，sylph effect 为 `0.283`（95% CI `0.140 to 0.427`），是该输出层级最强的被测工具。Chimera audited strict-read effect 为 `0.336`，数值上更大且更接近 published effect `0.384`。两者对应不同的输出对象，因此分别作为 sample-level abundance signal 和 audited read-level evidence 报告。

### Reference and Direct-alignment Robustness

Reference sensitivity 使用预先固定的 60-sample subset；样本选择不读取任何工具结果：每个队列从 paper C2-positive 样本的 abundance range 上等距选择 10 个，再按 disease/control 构成匹配 10 个 zero-Fna 样本，共 30 positive 和 30 zero-Fna。

审计参考库只在同一 clade 内删除 exact duplicate sequences，并保留跨 clade 相同序列以表达真实歧义。
原始 1,196 条 reference records 中删除 63 条同 clade exact duplicates，保留 1,133 条。
每个 read 分别取得每个 clade 内的最佳 alignment，再跨 clade 比较同一 `(nmatch, alignment length)` score；这一设计避免 combined-reference secondary-hit cap 决定哪些 clades 可见。

Primary robustness endpoint 仍为 all-positive `strict-C2-best reads/M AUC`：

| Method | AUC | Bootstrap 95% CI |
| --- | ---: | ---: |
| Chimera candidate audit, original combined panel | 0.867 | 0.767 to 0.942 |
| Chimera candidate audit, exact-dedup per-clade panel | **0.867** | **0.771 to 0.943** |
| Centrifuger candidate audit, exact-dedup per-clade panel | 0.751 | 0.621 to 0.863 |
| Kraken2 LF01 candidate audit, exact-dedup per-clade panel | 0.626 | 0.480 to 0.760 |
| Direct minimap2 alignment of all reads, exact-dedup per-clade panel | 0.792 | 0.671 to 0.896 |

在 `floor >=500` 时，对应 AUC 为 Chimera `0.994`、direct all-read alignment `0.949`、Centrifuger `0.894`、Kraken2 LF01 `0.697`。
Exact dedup 与 per-clade scoring 完整保留 Chimera 的结果；Chimera strict reads/M 的 median absolute change 仅 `2.17 reads/M`。
Direct all-read alignment 证明高信号样本本身可以被 targeted alignment 恢复，但在覆盖 published abundance range 的 60-sample subset 上，Chimera candidate evidence 给出了更清楚的 positive/zero-Fna separation。

### Reproducing the Statistical Tables

仓库保留了 760 个样本的标准化 manifest、四工具 sample-level signals 和三工具 read-audit metrics。以下命令只读取这些冻结输入，可重建主结果表；固定的 bootstrap 次数与随机种子会逐字段复现已提交的 TSV。

```bash
BASE=results/real/fna_c2_crc3_head3m
OUT=tmp/fna_c2_reproduction
mkdir -p "$OUT"

python scripts/fna_c2/make_detection_floor_tables.py \
  --sample-manifest "$BASE/sample_manifest.tsv" \
  --signal-table "$BASE/sample_level_signals.tsv" \
  --overall-out "$OUT/sample_level_overall_metrics.tsv" \
  --floor-out "$OUT/sample_level_detection_floor_auc.tsv"

python scripts/fna_c2/bootstrap_read_audit_auc.py \
  --sample-manifest "$BASE/sample_manifest.tsv" \
  --audit-table "$BASE/read_audit_sample_metrics.tsv" \
  --out "$OUT/tool_read_identity_audit/tool_audited_signal_auc_comparison.3tools.tsv" \
  --bootstrap 10000 --seed 20260710

python scripts/fna_c2/analyze_c2_specificity.py \
  --sample-manifest "$BASE/sample_manifest.tsv" \
  --audit-table "$BASE/read_audit_sample_metrics.tsv" \
  --signal-table "$BASE/sample_level_signals.tsv" \
  --out-dir "$OUT/specificity"

python scripts/fna_c2/analyze_crc_association.py \
  --metadata-table "$BASE/sample_manifest.tsv" \
  --combined-sample-table "$BASE/read_audit_sample_metrics.tsv" \
  --out-dir "$OUT/crc_read" --bootstrap 10000 --seed 20260710

python scripts/fna_c2/analyze_crc_sample_signal.py \
  --metadata-table "$BASE/sample_manifest.tsv" \
  --signal-table "$BASE/sample_level_signals.tsv" \
  --out-dir "$OUT/crc_signal" --bootstrap 10000 --seed 20260710
```

从原始 FASTQ 重新执行候选 read 提取、逐 clade alignment audit 和 60-sample robustness analysis 时，使用 `scripts/fna_c2/` 中相应的 audit、reference-build、sample-selection 和 sensitivity scripts。

### Result Files

| Path | Content |
| --- | --- |
| `fna_c2_crc3_head3m/sample_manifest.tsv` | 760-sample cohort, phenotype, accession, paper abundance, covariate and input-depth manifest |
| `fna_c2_crc3_head3m/sample_level_signals.tsv` | 760 samples x four tools, long-format native C2 signals for plotting |
| `fna_c2_crc3_head3m/read_audit_sample_metrics.tsv` | 760 samples x three read classifiers, candidate fate counts, rates and reads/M |
| `fna_c2_crc3_head3m/reference_manifest.tsv` | frozen reference input accession, clade, sequence count, bases and SHA256 |
| `fna_c2_crc3_head3m/sample_level_overall_metrics.tsv` | three-cohort overall sample-level C2 signal metrics |
| `fna_c2_crc3_head3m/sample_level_detection_floor_auc.tsv` | three-cohort sample-level detection-floor AUC |
| `fna_c2_crc3_head3m/tool_read_identity_audit/tool_audited_signal_auc_comparison.3tools.tsv` | three-tool read-level audit AUC, bootstrap CI and paired AUC differences |
| `fna_c2_crc3_head3m/c2_specificity_partial_association.tsv` | strict C2 evidence partial association with published C2 and C1 abundance |
| `fna_c2_crc3_head3m/crc_association_cohort.tsv` | cohort-level CRC association, bootstrap AUC and adjusted effects |
| `fna_c2_crc3_head3m/crc_association_meta.tsv` | three-cohort random-effects meta-analysis |
| `fna_c2_crc3_head3m/crc_sample_signal_association_cohort.tsv` | cohort-level CRC association for all four tools' native C2 signals |
| `fna_c2_crc3_head3m/crc_sample_signal_association_meta.tsv` | three-cohort native-signal meta-analysis including sylph |
| `fna_c2_crc3_head3m/robustness_sample_manifest.tsv` | pre-specified 60-sample reference-sensitivity subset |
| `fna_c2_crc3_head3m/robustness_sensitivity_auc.tsv` | original audit, exact-dedup per-clade audit and all-read direct-alignment AUC |
| `fna_c2_crc3_head3m/reference_panel_exact_dedup_summary.tsv` | per-clade exact-dedup accounting |

### Summary

Across three CRC cohorts, Chimera converts the published Fna C2 signal into auditable read-level evidence.
Across all C2-positive samples and at the detectable C2-signal floors, Chimera produces the strongest audited read-level separation among the read classifiers, with paired-bootstrap differences over Centrifuger and Kraken2.
The strict evidence tracks published C2 abundance after controlling for C1 abundance and recovers the CRC association in all three cohorts under the source paper's adjusted model.
The read-level result remains strongest after exact-reference deduplication, per-clade scoring and an all-read direct-alignment baseline.
