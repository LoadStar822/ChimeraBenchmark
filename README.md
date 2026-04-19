# ChimeraBenchmark

ChimeraBenchmark 汇集了论文中用于评估 Chimera 的宏基因组分类基准测试（metagenomic taxonomic benchmark）。它包含数据集说明、参考库规模、评估任务、结果表和复现脚本。

如果你是从论文来到这里，建议先阅读 `results/README.md`，再查看具体的 build、classify 和 profile 结果表。

## 目录

- [Benchmark 设计](#benchmark-设计)
- [数据与结果](#数据与结果)
- [复现方式](#复现方式)
- [仓库结构](#仓库结构)
- [方法边界](#方法边界)

## Benchmark 设计

本 benchmark 覆盖三类任务：

- 数据库构建（build）：比较不同工具在同一参考库上的构建成本。
- 逐读段分类（per-read classify）：比较每条 read 或 contig 在目标层级上的分类准确性。
- 丰度画像（profile）：比较样本层面的物种组成估计。

公开评估默认报告 `species` 和 `genus` 两个层级。逐读段结果使用真值映射，丰度结果使用样本组成真值。

## 数据与结果

- `results/README.md`：数据集、参考库、评估口径和指标说明。
- `results/builds/README.md`：数据库构建结果。
- `results/classify/README.md`：逐读段分类结果。
- `results/profile/README.md`：丰度画像结果。

`results/README.md` 中的两张 summary 表可直接对应论文补充材料：

- Reference Database Summary：参考库规模和来源。
- Benchmark Dataset Summary：实际输入数据、测序类型、读长和真值类型。

## 复现方式

查看可用命令：

```bash
python -m chimera_bench.cli -h
```

重新生成数据集与参考库 summary 表：

```bash
python -m chimera_bench.cli catalog --config configs --results-root results --resources-root resources
```

运行核心测试：

```bash
pytest -q tests/test_catalog.py tests/test_results_readme.py tests/test_cli.py
```

## 仓库结构

- `configs/`：数据集、参考库和工具运行配置。
- `chimera_bench/`：运行、汇总和评估代码。
- `results/`：结果表和各工具运行输出。
- `resources/`：统一分类体系、输入清单和 summary 表。
- `tests/`：README 生成、指标汇总和命令行行为测试。

## 方法边界

- 除线程数外，公开比较尽量使用各工具默认推断参数。
- 多步工具按端到端流程（end-to-end pipeline）统计运行时间和内存，例如 `kraken2 -> bracken`。
- 公开评估统一使用 NCBI 分类体系（NCBI taxonomy）快照 `resources/taxonomy/ncbi_20260408/`。
- 逐读段 F1 需要结合 `Truth Mapped Rate` 解读；当映射率不足时，结论只代表可映射子集。
