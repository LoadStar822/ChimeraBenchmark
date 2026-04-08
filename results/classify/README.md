# Classify Results

Auto-generated. Do not edit.

## Dataset: cami2-marine-long

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数合同的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 1327.807372 | 142.294285 | 0.826114 | 0.932455 | 0.868872 | 0.862185 | 0.865515 | 0.951177 | 0.967606 | 0.961332 | 0.957504 | 0.959414 |
| chimera | cami_refseq | 819.923478 | 79.227821 | 0.826114 | 0.958518 | 0.914434 | 0.88562 | 0.899796 | 0.951177 | 0.947996 | 0.976494 | 0.946621 | 0.961325 |
| ganon2 | cami_refseq | 4373.419369 | 212.943272 | 0.826114 | 0.785463 | 0.926892 | 0.783362 | 0.849105 | 0.951177 | 0.78398 | 0.986588 | 0.796933 | 0.881677 |
| kraken2 | cami_refseq | 243.954304 | 123.355232 | 0.826114 | 0.939645 | 0.867414 | 0.860703 | 0.864045 | 0.951177 | 0.976069 | 0.972072 | 0.968627 | 0.970346 |
| taxor | cami_refseq | 982.843509 | 86.356911 | 0.826114 | 0.794467 | 0.459832 | 0.388411 | 0.421115 | 0.951177 | 0.792985 | 0.470614 | 0.383871 | 0.42284 |

## Dataset: cami2-marine-long-sample0

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数合同的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 262.735984 | 142.185009 | 0.812604 | 0.928087 | 0.863819 | 0.856165 | 0.859975 | 0.931514 | 0.967862 | 0.963678 | 0.959891 | 0.961781 |
| ganon2 | cami_refseq | 967.721634 | 210.739708 | 0.812604 | 0.784991 | 0.93432 | 0.799042 | 0.861402 | 0.931514 | 0.782683 | 0.994387 | 0.810405 | 0.893018 |
| kraken2 | cami_refseq | 111.661773 | 122.170258 | 0.812604 | 0.935242 | 0.85727 | 0.849626 | 0.853431 | 0.931514 | 0.975731 | 0.97294 | 0.969533 | 0.971234 |
| taxor | cami_refseq | 109.52994 | 86.202934 | 0.812604 | 0.792311 | 0.442863 | 0.376194 | 0.406815 | 0.931514 | 0.790004 | 0.442864 | 0.36323 | 0.399113 |

## Dataset: cami2-marine-short

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数合同的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 1753.243818 | 141.592957 | 0.829501 | 0.839356 | 0.829513 | 0.819159 | 0.824304 | 0.959836 | 0.947514 | 0.955736 | 0.938697 | 0.94714 |
| chimera | cami_refseq | 718.197562 | 106.462322 | 0.829501 | 0.929008 | 0.916569 | 0.87813 | 0.896938 | 0.959836 | 0.918464 | 0.974812 | 0.91343 | 0.943123 |
| ganon2 | cami_refseq | 3722.749219 | 210.564232 | 0.829501 | 0.823832 | 0.921861 | 0.808752 | 0.86161 | 0.959836 | 0.822279 | 0.986277 | 0.82835 | 0.900441 |
| kraken2 | cami_refseq | 259.510867 | 123.464989 | 0.829501 | 0.833778 | 0.817576 | 0.806928 | 0.812217 | 0.959836 | 0.953209 | 0.963874 | 0.94673 | 0.955226 |
| taxor | cami_refseq | 8511.758909 | 107.677494 | 0.829501 | 0.829346 | 0.445245 | 0.389352 | 0.415427 | 0.959836 | 0.827794 | 0.447227 | 0.377656 | 0.409508 |

