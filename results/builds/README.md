# Database Build Results

> This page is maintained by hand and is not rewritten by new experiments; the machine-readable build table is `summary.tsv` in the same directory.

This page summarises the CAMI RefSeq builds used in the paper and the scaling experiment from Archaea to the full RefSeq.
The values used for the figures, together with tool versions, input sizes and run status, are recorded in `refseq_scaling.tsv`.

## Measurement

- The scaling experiment uses 32 threads throughout, on the same server with 2 x AMD EPYC 7763 and 1 TiB RAM.
- Apart from thread count, the frozen NCBI taxonomy path and output paths, all tools use their default build parameters; no fast-build options are used.
- `Elapsed` is the wall-clock time of all configured steps; `Peak RSS` is the peak resident memory reported by `/usr/bin/time -v`.
- `DB size` is the logical size of the final database files after temporary files are removed, excluding the `/usr/bin/time` resource logs.
- The default Chimera 1.8.0 database contains the core filter, the PROVE data, the presence sketch and the taxonomy.
- The sylph result consists of the default genome sketches without an external taxonomy file; it provides the sketches needed for sample-level profiling, not an index for per-read classification.

## CAMI RefSeq

This table lists the builds used in the paper; its machine-readable version is `summary.tsv` in the same directory.
It keeps the raw measurement fields. The tool versions are Centrifuger 1.1.0-r291, Chimera 1.8.0, ganon 2.1.0, Kraken2 2.1.3, Metabuli 1.2.0 and sylph 0.8.1.
The original `meta.json` of the sylph/CAMI RefSeq build was not kept, so its measurement is carried over from the earlier build table; the other five rows use the original run metadata.
The Chimera and ganon rows were measured by rebuilding with the same build configuration after the benchmark classification runs had finished, on a host shared with other jobs. The rebuilt databases were checked against the databases used for classification: Chimera's records are identical, with some files differing only in write order; ganon's input list and taxonomy files are byte-identical, and its index differs only in the input file paths it records. ganon builds a species-level database; see the tool notes in `results/README.md`.

| Tool | DB Name | Elapsed Seconds | Max RSS (KB) | DB Size | Started At | Finished At |
| --- | --- | ---: | ---: | ---: | --- | --- |
| centrifuger | cami_refseq | 937693.595364 | 675187264 | 141 GiB | 2026-01-20T12:32:50+08:00 | 2026-01-31T09:01:03+08:00 |
| chimera | cami_refseq | 3473.4634470939636 | 85550912 | 106 GiB | 2026-10-01T09:36:50+08:00 | 2026-10-01T10:34:43+08:00 |
| ganon | cami_refseq | 10845.479605674744 | 109470164 | 76.1 GiB | 2026-10-01T12:16:36+08:00 | 2026-10-01T15:17:22+08:00 |
| kraken2 | cami_refseq | 119277.976333 | 128058564 | 315 GiB | 2026-01-16T19:27:44+08:00 | 2026-01-18T04:35:42+08:00 |
| metabuli | cami_refseq | 34318.325497 | 87889432 | 242 GiB | 2026-08-21T00:52:46+08:00 | 2026-08-21T10:24:44+08:00 |
| sylph | cami_refseq | 1810.144336 | 27915272 | 21.6 GiB | 2026-01-02T12:49:53+08:00 | 2026-01-02T13:20:04+08:00 |

## RefSeq Scaling

The four datasets increase in reference size. Each row keeps the input snapshot, tool version and run status; exact values, parameter policies and target SHA-256 checksums are recorded in `refseq_scaling.tsv`.

Metabuli builds from an assembly list with `--gtdb 1`. To match the frozen `ncbi_20260408` taxonomy, only deleted taxids were remapped according to the current NCBI assembly metadata; the sequence files and the number of input entries are unchanged.

| Dataset | Snapshot | Input entries | Input size | Tool | Version | Build time | Peak RSS (GiB) | DB size (GiB) | Status |
| --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | --- |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | chimera | 1.8.0 | 10.701 s | 1.57 | 1.32 | measured |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | centrifuger | v1.0.12 | 4 min 07 s | 11.03 | 0.64 | measured |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | metabuli | 1.2.0 | 3 min 08 s | 70.16 | 4.52 | measured |
| Archaea | 2024-10-10 | 600 | 1.72 Gbp | sylph | 0.8.1 | 1.323 s | 0.45 | 0.06 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | kraken2 | 2.1.3 | 6 min 30 s | 26.70 | 2.10 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | bracken | 2.9 | 5 min 12 s | 3.00 | 2.50 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | ganon | 2.1.0 (IBF mode) | 23 s | 2.90 | 0.60 | measured |
| Archaea | 2024-10-10 | 601 | 1.73 Gbp | ganon2 | 2.1.0 | 53 s | 3.40 | 1.40 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | chimera | 1.8.0 | 4 min 53 s | 16.17 | 25.57 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | centrifuger | v1.0.12 | 4 h 39 min | 66.80 | 21.03 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | metabuli | 1.2.0 | 1 h 11 min | 91.87 | 115.07 | measured |
| CompleteONE | 2024-10-07 | 25,129 | 53.78 Gbp | sylph | 0.8.1 | 36.710 s | 3.33 | 1.95 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | kraken2 | 2.1.3 | 4 h 11 min | 66.60 | 62.70 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | bracken | 2.9 | 6 h 41 min | 66.20 | 99.40 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | ganon | 2.1.0 (IBF mode) | 18 min 55 s | 40.40 | 36.10 | measured |
| CompleteONE | 2024-09-26 | 25,091 | 52.32 Gbp | ganon2 | 2.1.0 | 38 min 07 s | 56.20 | 37.80 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | chimera | 1.8.0 | 19 min 43 s | 35.75 | 52.51 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | centrifuger | v1.0.12 | 34 h 22 min | 213.68 | 52.79 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | metabuli | 1.2.0 | 3 h 22 min | 87.06 | 135.87 | measured |
| Complete | 2024-10-07 | 58,671 | 189.36 Gbp | sylph | 0.8.1 | 2 min 43 s | 9.19 | 6.84 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | kraken2 | 2.1.3 | 14 h 31 min | 74.20 | 73.40 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | bracken | 2.9 | 9 h 15 min | 76.60 | 134.10 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | ganon | 2.1.0 (IBF mode) | 2 h 50 min | 991.10 | 322.20 | measured |
| Complete | 2024-10-07 | 58,672 | 189.36 Gbp | ganon2 | 2.1.0 | 2 h 05 min | 65.40 | 46.10 | measured |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | chimera | 1.8.0 | 8 h 26 min | 84.75 | 1,119.84 | measured |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | centrifuger | v1.0.12 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | metabuli | 1.2.0 | 71 h 44 min | 403.31 | - | OOM |
| RefSeq | 2024-10-16 | 403,263 | 3.43 Tbp | sylph | 0.8.1 | 42 min 00 s | 130.25 | 108.69 | measured |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | kraken2 | 2.1.3 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | bracken | 2.9 | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | ganon | 2.1.0 (IBF mode) | - | - | - | OOM |
| RefSeq | 2024-10-16 | 403,253 | 3.43 Tbp | ganon2 | 2.1.0 | - | - | - | OOM |

### Scaling Result

- On the three levels that the read classifiers completed, Chimera builds 23.1x, 57.1x and 104.6x faster than Centrifuger.
- At the Complete level, Chimera's peak RSS is 35.75 GiB against 213.68 GiB for Centrifuger, about 1/6.0 of it.
- Chimera completes the 3.43 Tbp full RefSeq with default parameters in 8 h 26 min at 84.75 GiB peak RSS. The final 1.09 TiB database consists of a 542.46 GiB core filter, 553.18 GiB of PROVE data (`local` in the `db_size_components` column), a 23.75 GiB presence sketch and the taxonomy.
- Centrifuger and Metabuli complete Archaea, CompleteONE and Complete; their RefSeq builds fail with OOM.
- sylph is faster and smaller at every shared level; its result corresponds to genome sketch construction. Chimera's build advantage concerns a full index that supports per-read classification and PROVE.
