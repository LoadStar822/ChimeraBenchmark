# ChimeraBenchmark (Legacy Snapshot)

This repository currently archives an older benchmarking implementation.

- All legacy code and documents live under `legacy/`.
- New benchmarking code has not been rebuilt yet.

If you are looking for the previous scripts, see `legacy/`.

## 实验公平性与默认参数策略

为保证论文实验公平性，本项目遵循以下原则：

- 统一使用 64 线程运行所有软件。
- 除非软件将 classify 与 profile 完全分离（必须额外生成 profile 输出），否则尽量使用软件默认参数。
- 仅为**记录评估结果**而增加输出开关（例如 per-read 输出），不改变算法行为与默认推断逻辑。

### Ganon（ganon2）当前 benchmark 的默认跑法

保持 ganon 的默认推断逻辑（含默认 EM 多重命中处理），仅增加输出以便评估：

- classify 仅追加 `--output-one`、`--output-unclassified`、`--skip-report`（不改变算法结果）。
- 不显式设置 `--multiple-matches / --rel-cutoff / --rel-filter / --fpr-query`，全部使用默认值。
- profile 由 `ganon report` 从 `.rep` 生成 `.tre`（reads 与 abundance 两类）。

如需变更该策略，必须在 README 中明确记录，并注明原因。
