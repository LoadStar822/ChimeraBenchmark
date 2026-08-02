# Database Build Results

本页汇总论文正式 CAMI RefSeq 构建结果，以及从 Archaea 到完整 RefSeq 的递增规模实验。
用于作图的数值、工具版本、输入规模和运行状态统一记录在 `refseq_scaling.tsv`。

## Measurement

- 递增规模实验统一使用 32 threads，运行于 2 x AMD EPYC 7763、1 TiB RAM 的同一台服务器。
- 同组工具除线程数、冻结的 NCBI taxonomy 路径和输出路径外均使用默认构建参数；没有使用 fast-build 参数。
- `Elapsed` 是完整配置步骤的 wall-clock time；`Peak RSS` 是 `/usr/bin/time -v` 记录的峰值常驻内存。
- `DB size` 是临时文件清理后最终数据库产物的逻辑大小，不包含 `/usr/bin/time` 资源日志。
- Chimera 1.7.1 的默认数据库包含 core filter、taxonomy 和 local read-resolution data。
- sylph 的结果是默认 genome sketches，不包含外部 taxonomy 文件；它提供 sample-level profiling 所需的 sketch，而不是逐 read 分类索引。

## CAMI RefSeq

这是 `results/builds/summary.tsv` 中进入论文矩阵的正式构建结果。
首表保留 benchmark freeze 所需的原始字段；对应工具版本为 Centrifuger 1.1.0-r291、Chimera 1.7.1、Kraken2 2.1.3 和 sylph 0.8.1。

| Tool | DB Name | Elapsed Seconds | Max RSS (KB) | DB Size | Started At | Finished At |
| --- | --- | ---: | ---: | ---: | --- | --- |
| centrifuger | cami_refseq | 937693.595364 | 675187264 | 141 GiB | 2026-01-20T12:32:50+08:00 | 2026-01-31T09:01:03+08:00 |
| chimera | cami_refseq | 3387.0917110443115 | 85313256 | 103 GiB | 2026-07-24T17:46:51+08:00 | 2026-07-24T18:43:18+08:00 |
| kraken2 | cami_refseq | 119277.976333 | 128058564 | 315 GiB | 2026-01-16T19:27:44+08:00 | 2026-01-18T04:35:42+08:00 |
| sylph | cami_refseq | 1810.144336 | 27915272 | 21.6 GiB | 2026-01-02T12:49:53+08:00 | 2026-01-02T13:20:04+08:00 |

## RefSeq Scaling

四级数据集按参考规模递增。每行保留对应的输入快照、工具版本和运行状态；精确数值、参数政策及 target SHA-256 记录在 `refseq_scaling.tsv`。

| Dataset | Snapshot | Input entries | Input size | Tool | Version | Build time | Peak RSS (GiB) | DB size (GiB) | Status |
| --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | --- |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | chimera | 1.7.1 (f888829) | 11.829 s | 1.57 | 1.29 | measured |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | centrifuger | v1.0.12-10-gdf5bd32 | 4 min 07 s | 11.03 | 0.64 | measured |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | sylph | 0.8.1 | 1.323 s | 0.45 | 0.06 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | kraken2 | 2.1.3 | 6 min 30 s | 26.70 | 2.10 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | bracken | 2.9 | 5 min 12 s | 3.00 | 2.50 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | ganon | 2.1.0 (IBF mode) | 23 s | 2.90 | 0.60 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | ganon2 | 2.1.0 | 53 s | 3.40 | 1.40 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | taxor | 0.1.3 | 1 min 19 s | 1.80 | 0.40 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | chimera | 1.7.1 (f888829) | 4 min 11 s | 16.17 | 24.98 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | centrifuger | v1.0.12-10-gdf5bd32 | 4 h 39 min | 66.80 | 21.03 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | sylph | 0.8.1 | 36.710 s | 3.33 | 1.95 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | kraken2 | 2.1.3 | 4 h 11 min | 66.60 | 62.70 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | bracken | 2.9 | 6 h 41 min | 66.20 | 99.40 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | ganon | 2.1.0 (IBF mode) | 18 min 55 s | 40.40 | 36.10 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | ganon2 | 2.1.0 | 38 min 07 s | 56.20 | 37.80 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | taxor | 0.1.3 | 3 h 56 min | 18.60 | 12.20 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | chimera | 1.7.1 (f888829) | 15 min 08 s | 35.75 | 51.35 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | centrifuger | v1.0.12-10-gdf5bd32 | 34 h 22 min | 213.68 | 52.79 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | sylph | 0.8.1 | 2 min 43 s | 9.19 | 6.84 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | kraken2 | 2.1.3 | 14 h 31 min | 74.20 | 73.40 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | bracken | 2.9 | 9 h 15 min | 76.60 | 134.10 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | ganon | 2.1.0 (IBF mode) | 2 h 50 min | 991.10 | 322.20 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | ganon2 | 2.1.0 | 2 h 05 min | 65.40 | 46.10 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | taxor | 0.1.3 | 6 h 30 min | 38.50 | 32.10 | measured |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | chimera | 1.7.1 (f888829) | 8 h 54 min | 84.76 | 1,096.08 | measured |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | centrifuger | v1.0.12-10-gdf5bd32 | - | - | - | not run |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | sylph | 0.8.1 | 42 min 00 s | 130.25 | 108.69 | measured |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | kraken2 | 2.1.3 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | bracken | 2.9 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | ganon | 2.1.0 (IBF mode) | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | ganon2 | 2.1.0 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | taxor | 0.1.3 | - | - | - | OOM |

### Scaling Result

- 在三个实际完成的 read-classifier 层级，Chimera 相对 Centrifuger 分别快 20.9x、66.6x 和 136.3x。
- Complete 层级中，Chimera peak RSS 为 35.75 GiB，Centrifuger 为 213.68 GiB；Chimera 约为后者的 1/6.0。
- Chimera 使用默认参数完成 3.43 Tbp Full RefSeq：8 h 54 min、84.76 GiB peak RSS。最终 1.07 TiB 数据库由 542.44 GiB core filter、553.18 GiB local read-resolution data 和 taxonomy 构成。
- sylph 在所有共有层级都更快且更小；该结果对应 genome sketch construction。Chimera 的构建优势定位于支持逐 read 分类和 local resolution 的完整索引。
