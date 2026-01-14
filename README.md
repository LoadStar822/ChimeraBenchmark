# ChimeraBenchmark

本仓库用于运行与汇总 Chimera 的基准测试（build / classify / profile）。

- 所有结果统一写入 `results/`（默认不提交原始输出与日志，仅提交聚合后的 README）。
- 实验与数据集配置在 `configs/` 下。
- 旧版脚本与文档已移动到 `legacy` 分支（见 `legacy` branch）。

## 实验公平性与默认参数策略

为保证论文实验公平性，本项目遵循以下原则：

- 统一使用相同线程数运行所有软件（默认 **32**，除非另行说明）。
- 除非软件将 classify 与 profile 完全分离（必须额外生成 profile 输出），否则尽量使用软件默认参数。
- 仅为**记录评估结果**而增加输出开关（例如 per-read 输出），不改变算法行为与默认推断逻辑。

### Ganon（ganon2）当前 benchmark 的默认跑法

保持 ganon 的默认推断逻辑（含默认 EM 多重命中处理），仅增加输出以便评估：

- classify 仅追加 `--output-one`、`--output-unclassified`、`--skip-report`（不改变算法结果）。
- 不显式设置 `--multiple-matches / --rel-cutoff / --rel-filter / --fpr-query`，全部使用默认值。
- profile 由 `ganon report` 从 `.rep` 生成 `.tre`（reads 与 abundance 两类）。

如需变更该策略，必须在 README 中明确记录，并注明原因。
