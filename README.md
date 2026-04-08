# ChimeraBenchmark

本仓库用于运行与汇总 Chimera 的基准测试（build / classify / profile）。

- 所有结果统一写入 `results/`（默认不提交原始输出与日志，仅提交聚合后的 README）。
- 实验与数据集配置在 `configs/` 下。
- 旧版脚本与文档已移动到 `legacy` 分支（见 `legacy` branch）。

## 实验公平性与默认参数策略

为保证论文实验公平性，本项目遵循以下原则：

- 统一使用相同线程数运行所有软件（默认 **32**，除非另行说明）。
- 推断参数保持软件默认；本仓库不再接受为了“好看结果”而改推断策略。
- 仅允许为了**导出评估结果**而增加输出开关；这类开关不能改变算法行为与默认推断逻辑。
- 若某工具必须依赖多步 pipeline（例如 `kraken2 -> bracken`），runtime / memory 必须按端到端 pipeline 解释。
- 评估时统一使用同一份 NCBI taxonomy 快照，不再混用各工具 build 目录里的 rank / name 定义。

### Ganon（ganon2）当前 benchmark 的默认跑法

保持 ganon 的默认推断逻辑，仅增加输出以便评估：

- classify 仅追加 `--output-one`、`--output-unclassified`、`--skip-report`（不改变算法结果）。
- 不显式设置 `--multiple-matches / --rel-cutoff / --rel-filter / --fpr-query`，全部使用默认值。
- profile 由 `ganon report` 从 `.rep` 生成 `.tre`（reads 与 abundance 两类）。

如需变更该策略，必须在 README 中明确记录，并注明原因。

## 历史结果刷新状态

- 当前仓库里已落盘的 `bracken` 结果，已按新的端到端 `kraken2 -> bracken` 合同刷新。
- 当前仓库里已落盘的 `taxor` 结果，包含旧的 benchmark 特殊处理痕迹，刷新后才能用于论文主表。
- 当前仓库里已落盘的 `ganon` 结果，已按新的默认参数合同刷新。

## 统一 taxonomy 快照

- 当前评估统一使用官方 NCBI taxonomy 快照：
  `results/taxonomy/ncbi_20260408/`
- 评估依赖文件：
  - `nodes.dmp`
  - `names.dmp`
  - `source.json`
  - `SHA256SUMS`
- 这份快照用于统一 rank 判定和名称归一化，避免不同工具各自 `.tax` 的差异直接进入评估口径。

## Per-read 主指标合同

- `species exact` 只统计能提升到 `species` 的真值，不要求原始真值的 rank 恰好等于 `species`。
- `genus exact` 只统计能提升到 `genus` 的真值，不要求原始真值的 rank 恰好等于 `genus`。
- 若原始真值是 `strain/subspecies/isolate`，但能唯一上提到某个 `species/genus`，则进入对应主指标分母。
- 若原始真值只有 `family/order/class`，无法上提到更细层级，则不进入该层级主指标分母。
- 公开结果必须同时展示 `Truth Mapped Rate`，否则 `species/genus F1` 不可单独引用。
