# Classify Results

Auto-generated. Do not edit.

## Dataset: cami-strain-madness-long

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 13575.142525 | 144.22858 | 0.908648 | 0.494245 | 0.493303 | 0.414794 | 0.450655 | 0.917533 | 0.638937 | 0.693488 | 0.582371 | 0.633091 |
| chimera | cami_refseq | 6242.094637 | 72.904514 | 0.908648 | 0.892819 | 0.065646 | 0.058551 | 0.061896 | 0.917533 | 0.622109 | 0.06763 | 0.060323 | 0.063768 |
| ganon2 | cami_refseq | 62372.79518 | 210.946983 | 0.908648 | 0 | 0.1 | 0 | 0 | 0.917533 | 0 | 0.2 | 0 | 0 |
| kraken2 | cami_refseq | 4312.531238 | 126.772358 | 0.908648 | 0.400711 | 0.382187 | 0.268151 | 0.315171 | 0.917533 | 0.506199 | 0.553196 | 0.387473 | 0.455736 |

## Dataset: cami-strain-madness-long-sample0

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 4614.12194 | 141.858158 | 0.888111 | 0.473686 | 0.473697 | 0.395524 | 0.431095 | 0.965477 | 0.643029 | 0.715039 | 0.583446 | 0.642574 |
| chimera | cami_refseq | 154.508462 | 72.282906 | 0.888111 | 0.895223 | 0.068762 | 0.061505 | 0.064932 | 0.965477 | 0.617912 | 0.067543 | 0.060472 | 0.063812 |
| ganon2 | cami_refseq | 5612.797757 | 210.805557 | 0.888111 | 0 | 0 | 0 | 0 | 0.965477 | 0 | 0 | 0 | 0 |
| kraken2 | cami_refseq | 399.296068 | 120.740246 | 0.888111 | 0.383283 | 0.3744 | 0.257906 | 0.305422 | 0.965477 | 0.506943 | 0.595465 | 0.398151 | 0.477216 |

## Dataset: cami-strain-madness-short

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。
集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数；`Elapsed` 为 sample 运行时间总和，`Max RSS` 为最大值，read 数为总和；率值和准确率为 sample 算术平均。

| Tool | DB | Samples | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq |  | 46035.597943 | 141.273022 | 0.921763 | 0.852283 | 0.825675 | 0.82428 | 0.824977 | 0.929726 | 0.901 | 0.893789 | 0.89193 | 0.892859 |
| chimera | cami_refseq |  | 56649.267445 | 72.209198 | 0.921763 | 0.994036 | 0.934087 | 0.931152 | 0.932617 | 0.929726 | 0.992255 | 0.953553 | 0.949119 | 0.951331 |
| ganon2 | cami_refseq | 50 | 433823.63865 | 210.488087 | 0.921762 | 0.107288 | 0.95924 | 0.103576 | 0.186961 | 0.929728 | 0.107288 | 0.989733 | 0.10647 | 0.192255 |
| kraken2 | cami_refseq |  | 3788.487527 | 130.382889 | 0.921763 | 0.614644 | 0.529935 | 0.529 | 0.529467 | 0.929726 | 0.787993 | 0.754908 | 0.753273 | 0.75409 |

## Dataset: cami-strain-madness-short-sample0

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 1552.178563 | 141.259415 | 0.922091 | 0.833328 | 0.795925 | 0.794741 | 0.795333 | 0.922742 | 0.878265 | 0.859314 | 0.858032 | 0.858672 |
| chimera | cami_refseq | 790.460696 | 71.976997 | 0.922091 | 0.996835 | 0.930317 | 0.927503 | 0.928908 | 0.922742 | 0.995993 | 0.945389 | 0.942494 | 0.94394 |
| ganon2 | cami_refseq | 4936.976575 | 210.249363 | 0.922091 | 0.107915 | 0.965146 | 0.104256 | 0.188183 | 0.922742 | 0.107915 | 0.98769 | 0.106665 | 0.192537 |
| kraken2 | cami_refseq | 135.153647 | 121.81501 | 0.922091 | 0.605807 | 0.518509 | 0.517703 | 0.518105 | 0.922742 | 0.762302 | 0.715065 | 0.71395 | 0.714507 |

## Dataset: cami2-marine-long

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 1327.807372 | 142.294285 | 0.826114 | 0.932455 | 0.868872 | 0.862185 | 0.865515 | 0.951177 | 0.967606 | 0.961332 | 0.957504 | 0.959414 |
| chimera | cami_refseq | 266.858838 | 73.034859 | 0.826114 | 0.978327 | 0.901204 | 0.888053 | 0.89458 | 0.951177 | 0.961255 | 0.964128 | 0.947341 | 0.955661 |
| ganon2 | cami_refseq | 4373.419369 | 212.943272 | 0.826114 | 0.785463 | 0.926892 | 0.783362 | 0.849105 | 0.951177 | 0.78398 | 0.986588 | 0.796933 | 0.881677 |
| kraken2 | cami_refseq | 243.954304 | 123.355232 | 0.826114 | 0.939645 | 0.867414 | 0.860703 | 0.864045 | 0.951177 | 0.976069 | 0.972072 | 0.968627 | 0.970346 |

## Dataset: cami2-marine-long-sample0

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 262.735984 | 142.185009 | 0.812604 | 0.928087 | 0.863819 | 0.856165 | 0.859975 | 0.931514 | 0.967862 | 0.963678 | 0.959891 | 0.961781 |
| ganon2 | cami_refseq | 967.721634 | 210.739708 | 0.812604 | 0.784991 | 0.93432 | 0.799042 | 0.861402 | 0.931514 | 0.782683 | 0.994387 | 0.810405 | 0.893018 |
| kraken2 | cami_refseq | 111.661773 | 122.170258 | 0.812604 | 0.935242 | 0.85727 | 0.849626 | 0.853431 | 0.931514 | 0.975731 | 0.97294 | 0.969533 | 0.971234 |

## Dataset: cami2-marine-short

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。

| Tool | DB | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 1753.243818 | 141.592957 | 0.829501 | 0.839356 | 0.829513 | 0.819159 | 0.824304 | 0.959836 | 0.947514 | 0.955736 | 0.938697 | 0.94714 |
| chimera | cami_refseq | 497.593706 | 72.679955 | 0.829501 | 0.963676 | 0.909019 | 0.892587 | 0.900728 | 0.959836 | 0.950959 | 0.964393 | 0.934503 | 0.949213 |
| ganon2 | cami_refseq | 3722.749219 | 210.564232 | 0.829501 | 0.823832 | 0.921861 | 0.808752 | 0.86161 | 0.959836 | 0.822279 | 0.986277 | 0.82835 | 0.900441 |
| kraken2 | cami_refseq | 259.510867 | 123.464989 | 0.829501 | 0.833778 | 0.817576 | 0.806928 | 0.812217 | 0.959836 | 0.953209 | 0.963874 | 0.94673 | 0.955226 |

## Dataset: prjna637878-supported19

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。
集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数；`Elapsed` 为 sample 运行时间总和，`Max RSS` 为最大值，read 数为总和；率值和准确率为 sample 算术平均。

| Tool | DB | Samples | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chimera | cami_refseq | 19 | 19716.048342 | 72.358036 | 1 | 0.9106 | 0.721678 | 0.664895 | 0.69019 | 0.992838 | 0.893599 | 0.805712 | 0.742209 | 0.770534 |
| ganon2 | cami_refseq | 16 | 94590.988172 | 294.623417 | 1 | 0 | 0 | 0 | 0 | 0.991495 | 0 | 0 | 0 | 0 |
| kraken2 | cami_refseq | 19 | 6344.123134 | 124.525604 | 1 | 0.639691 | 0.499353 | 0.472261 | 0.484409 | 0.992645 | 0.840341 | 0.735709 | 0.692797 | 0.711769 |

## Dataset: prjna637878-supported19-single-read

### Per-read Metrics

本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。
例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。
原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。
Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。
只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。
集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数；`Elapsed` 为 sample 运行时间总和，`Max RSS` 为最大值，read 数为总和；率值和准确率为 sample 算术平均。

| Tool | DB | Samples | Elapsed (s) | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) | Truth Mapped Rate (genus) | Pred Mapped Rate (genus) | Precision (genus) | Recall (genus) | F1 (genus) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 19 | 30255.753284 | 141.351318 | 1 | 0.598278 | 0.465997 | 0.43533 | 0.448398 | 0.992453 | 0.823077 | 0.731804 | 0.681198 | 0.702425 |
| chimera | cami_refseq | 19 | 24227.181703 | 73.901546 | 1 | 0.901788 | 0.708699 | 0.644661 | 0.672324 | 0.992453 | 0.885961 | 0.799311 | 0.726522 | 0.75803 |
| kraken2 | cami_refseq | 19 | 6588.018495 | 123.828533 | 1 | 0.576879 | 0.451795 | 0.419893 | 0.43352 | 0.992453 | 0.811111 | 0.723497 | 0.670724 | 0.692936 |

