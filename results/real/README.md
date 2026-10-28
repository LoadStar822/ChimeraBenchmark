# Real-data Results

This directory holds the real-cohort experiments outside the default benchmark.
The current experiment tests, on three published CRC cohorts, whether tools can recover the Fna C2 signal (that is, *F. animalis*) reported in the original study (Zepeda-Rivera et al., Nature 628, 424–432, 2024) with a species-level *Fusobacterium* reference panel, and turns the sample-level signal into auditable read-level evidence.
Results are reported separately for the sample-level signal from each tool's native output, the read-level audit and the CRC association, matching the native output level of each tool.

> This page describes the experiment design and interprets the results. It is maintained by hand and is not rewritten by new experiments.
> The result tables are in `results/real/fna_c2_crc3_head3m/`; the file list is at the end of this page.

## Dataset: fna-c2-crc3-head3m

### Experiment Design

The Fna C2 CRC stool metagenome experiment combines the Yachida, Yu and Wirbel cohorts, 760 R1 samples in total (392 CRC, 368 control).
Each sample uses the first 3M R1 reads; 2 samples whose source files have fewer than 3M reads use all available R1 reads.

Sample labels come from the Fna C1/C2 abundance signal reported in the original study and are not redefined by Chimera or any other tool:

| Role | Samples | Meaning | Included in C2-positive vs zero-Fna AUC |
| --- | ---: | --- | --- |
| paper C2-positive | 137 | Fna C2 signal reported in the original study | yes |
| paper zero-Fna | 619 | zero Fna signal reported in the original study | yes |
| paper Fna-positive / C2-zero | 4 | Fna signal but zero C2 reported in the original study | no |

### Reference Panel

The classification reference is a frozen species-level *Fusobacterium* panel: 188 reference genomes from 8 species, 1,161 sequences and 451,155,567 bp in total (`reference_manifest.tsv`).
Of these, 147 genomes come from the three clade genome files of Supplementary Table 1 of the original study, split into single-strain genomes; the other 41 are closely related NCBI genomes (40 RefSeq, 1 GenBank).
Only one copy of each strain (same BioSample or same strain name) is kept; the 24 removed genomes are listed in `species_panel_dropped.tsv`.

The species of each genome was fixed by a set rule before any classification run: when the NCBI ANI taxonomy check reports species_match, derived_species_match or mismatch and the species of the closest type strain exists in the NCBI taxonomy snapshot (2026-04-08), that species is used (171 genomes); otherwise the NCBI organism name is used (16); when neither is available, the subspecies annotation of the original study is used (1).
Under this rule, all 65 genomes labelled Fna C2 in the original study fall into *F. animalis* and all 24 labelled Fna C1 into *F. watanabei*, in agreement with the NCBI ANI check and with the species delineation of Sivertsen et al. (mBio 2025, doi:10.1128/mbio.00941-25).

| Species | NCBI taxid | Signal category | Genomes |
| --- | ---: | --- | ---: |
| *F. animalis* | 76859 | Fna C2 | 67 |
| *F. watanabei* | 2686067 | Fna C1 | 33 |
| *F. polymorphum* | 76857 | non-C1/C2 | 40 |
| *F. nucleatum* | 851 | non-C1/C2 | 18 |
| *F. vincentii* | 155615 | non-C1/C2 | 13 |
| *F. mortiferum* | 850 | other | 8 |
| *F. varium* | 856 | other | 8 |
| *F. paranimalis* | 3453848 | other | 1 |

A read that a tool assigns to *F. animalis* counts as C2 signal, and one assigned to *F. watanabei* counts as C1 signal.

The read-level audit uses a separate frozen paper-clade audit reference (`audit_reference_manifest.tsv`): 68 reference files with 1,196 sequences and 507,495,341 bp, comprising the three clade files of the original study and 65 NCBI genomes (64 RefSeq, 1 GenBank), labelled with the clades of the original study.
The audit reference was frozen before the species-level panel and does not change with the classification panel.

### Tool Settings

All three cohorts use the same fixed input reads, and the six tools build their databases from the same species-level panel. Apart from thread count, output paths and the settings noted below, the build and the classification or profiling steps use each tool's default parameters.

| Tool | Version | Run setting |
| --- | --- | --- |
| Chimera | 1.8.0 | default `build` and `classify`, 8 threads per sample; `--profile-cami` only requests an additional output and does not change read classification |
| Centrifuger | source `v1.0.12-10-gdf5bd32` | default single-read classification, 8 threads |
| Kraken2 | 2.1.3 | default classification, 8 threads; the DB build uses load factor 0.1, the same build setting as for the earlier panel (on this panel both the default and the 0.2 builds hit a compact hash table capacity failure) |
| Metabuli | 1.2.0 | default single-end classification (`--seq-mode 1`; the input is R1 only), 8 threads |
| ganon | 2.1.0 | species-level database (`build-custom --level species`, the default level of `ganon build`; `--skip-genome-size` only affects abundance reports); per-read classification uses the `--binning` preset that the ganon documentation gives for read binning, 8 threads |
| sylph | 0.8.1 | default genome sketch parameter `c=200` and default `profile`, 8 threads |
| minimap2 audit | 2.24-r1122 | orthogonal short-read alignment audit (`-x sr`) |

### Metric Definitions

`expected C2 reads` converts the Fna C2 abundance of the original study to the input depth of this experiment:

```text
expected C2 reads = paper Fna C2 percent * input reads / 100
```

A `detection floor` keeps as positives only the paper C2-positive samples whose `expected C2 reads` reach the given floor; zero-Fna samples are unchanged.
For example, `floor >=1000` evaluates only positive samples in which, at the current sampling depth, the signal of the original study is expected to contribute at least 1000 C2 reads.

The `sample-level C2 signal` is defined by the type of tool output:

| Tool type | Signal definition |
| --- | --- |
| Chimera / Centrifuger / Kraken2 LF01 / Metabuli / ganon | reads assigned to *F. animalis* (Fna C2) per million input reads |
| sylph | sequence abundance percentage summed over *F. animalis* references |

`AUC` measures whether a tool's signal ranks paper C2-positive samples above zero-Fna samples.
Higher values mean a clearer separation of positive and negative samples; `0.5` is close to random ranking and `1.0` is complete separation.

`Pearson vs paper C2` is the linear correlation between a tool's C2 signal and the Fna C2 abundance signal of the original study.

The `read-level audit` checks the C2 candidate reads called by a tool with an orthogonal read alignment.
Candidate reads (reads a tool assigns to *F. animalis*) are realigned to the audit reference and classified by their best-scoring alignments:

| Audit metric | Definition |
| --- | --- |
| raw candidate reads/M | reads per million that the tool calls Fna C2, i.e. the native C2 signal of a read-level tool |
| audited C2-supported reads/M | candidate reads per million whose best-scoring alignments include an Fna C2 reference, including strict and tied support |
| audited strict-C2-best reads/M | candidate reads per million whose best-scoring alignments fall only on Fna C2 references |
| call precision (strict-C2-best rate) | the fraction of candidate reads that are strict-C2-best |

The primary read-level endpoint is `audited strict-C2-best reads/M`. The tables below summarise the main results of the combined three-cohort analysis.

### Sample-level Signals

The following results compute the sample-level C2 signal from each tool's native output.

| Tool | Signal unit | Samples | C2-positive mean | zero-Fna mean | all-positive AUC | Pearson vs paper C2 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | reads/M | 760 | 533.156 | 30.749 | **0.865** | **0.951** |
| ganon | reads/M | 760 | 865.119 | 340.556 | 0.761 | 0.907 |
| Centrifuger | reads/M | 760 | 716.725 | 352.998 | 0.727 | 0.923 |
| Kraken2 LF01 | reads/M | 760 | 186.895 | 72.342 | 0.683 | 0.883 |
| Metabuli | reads/M | 760 | 1374.961 | 1074.412 | 0.627 | 0.833 |
| sylph | sequence abundance % | 760 | 17.498 | 0.000 | 0.624 | 0.565 |

| Tool | all positives (n=137) | floor >=500 (n=50) | floor >=1000 (n=38) | floor >=2000 (n=26) | floor >=5000 (n=12) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.865 | 0.986 | 0.993 | 0.995 | 0.993 |
| ganon | 0.761 | 0.920 | 0.968 | 0.988 | 0.992 |
| Centrifuger | 0.727 | 0.898 | 0.957 | 0.994 | 0.996 |
| Kraken2 LF01 | 0.683 | 0.877 | 0.932 | 0.987 | 0.981 |
| Metabuli | 0.627 | 0.790 | 0.874 | 0.971 | 0.986 |
| sylph | 0.624 | 0.840 | 0.908 | 0.962 | 0.958 |

Chimera's native C2 signal reaches an AUC of `0.865` over all positive samples, `0.104` above the next best, ganon (`0.761`; paired bootstrap 95% CI `0.058 to 0.149`), and has the highest Pearson correlation with the C2 abundance of the original study (`0.951`).
Among read-level tools, Chimera has the lowest background in zero-Fna samples: mean `30.7` reads/M and median `4.7` reads/M; the means of Centrifuger, Kraken2 LF01, Metabuli and ganon are `353.0`, `72.3`, `1074.4` and `340.6` reads/M.

### Read-level Audit

The read-level audit projects the C2 candidate reads called by each tool onto alignment-supported C2 evidence.
The same AUC compares positive samples with zero-Fna samples, with the signal changed from raw calls to audited reads/M.
`raw candidate reads/M` equals the native signal of the read-level tools in the previous section and is not repeated.

#### Audited Strict-C2-best Reads per Million

| Tool | all positives | floor >=500 | floor >=1000 | floor >=2000 | floor >=5000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chimera | 0.897 | 0.996 | 0.999 | 0.999 | 0.999 |
| Centrifuger | 0.750 | 0.928 | 0.971 | 0.997 | 0.999 |
| Kraken2 LF01 | 0.642 | 0.824 | 0.891 | 0.966 | 0.971 |
| Metabuli | 0.896 | 0.990 | 0.996 | 0.998 | 0.999 |
| ganon | 0.723 | 0.931 | 0.972 | 0.997 | 0.998 |

A paired bootstrap with 10,000 replicates, stratified by cohort and sample role, gives the confidence intervals below. `Chimera minus tool` is computed on the same resampled samples in each replicate and therefore tests the AUC difference on the same samples directly.

| Positive set | Tool | AUC (95% CI) | Chimera minus tool (95% CI) | Two-sided P |
| --- | --- | ---: | ---: | ---: |
| all positives | Chimera | **0.897 (0.862 to 0.929)** | - | - |
| all positives | Centrifuger | 0.750 (0.703 to 0.795) | **0.148 (0.102 to 0.193)** | < 2e-4 |
| all positives | Kraken2 LF01 | 0.642 (0.592 to 0.692) | **0.255 (0.200 to 0.308)** | < 2e-4 |
| all positives | Metabuli | 0.896 (0.867 to 0.922) | 0.001 (-0.029 to 0.028) | 0.896 |
| all positives | ganon | 0.723 (0.675 to 0.772) | **0.174 (0.127 to 0.222)** | < 2e-4 |
| floor >=500 | Chimera | **0.996 (0.992 to 0.999)** | - | - |
| floor >=500 | Centrifuger | 0.928 (0.887 to 0.964) | **0.068 (0.033 to 0.108)** | < 2e-4 |
| floor >=500 | Kraken2 LF01 | 0.824 (0.758 to 0.887) | **0.172 (0.110 to 0.237)** | < 2e-4 |
| floor >=500 | Metabuli | 0.990 (0.984 to 0.996) | **0.0057 (0.0021 to 0.0103)** | 0.0010 |
| floor >=500 | ganon | 0.931 (0.890 to 0.967) | **0.065 (0.030 to 0.106)** | < 2e-4 |
| floor >=1000 | Chimera | **0.999 (0.996 to 1.000)** | - | - |
| floor >=1000 | Centrifuger | 0.971 (0.941 to 0.996) | **0.028 (0.003 to 0.057)** | 0.0010 |
| floor >=1000 | Kraken2 LF01 | 0.891 (0.828 to 0.947) | **0.108 (0.052 to 0.170)** | < 2e-4 |
| floor >=1000 | Metabuli | 0.996 (0.992 to 0.999) | **0.0026 (0.0004 to 0.0054)** | 0.0138 |
| floor >=1000 | ganon | 0.972 (0.944 to 0.994) | **0.027 (0.005 to 0.054)** | < 2e-4 |

Over all positive samples, Chimera's strict-C2-best AUC is clearly higher than those of Centrifuger, Kraken2 LF01 and ganon, and level with the audited Metabuli; at `floor >=500` and `floor >=1000`, Chimera is higher than all four tools, and none of the 95% CIs of the paired differences contains 0.
Metabuli's native signal reaches an AUC of only `0.627` and rises to `0.896` only after the orthogonal alignment removes most background calls; Chimera's native signal reaches `0.865` without any audit.
Confidence intervals for the secondary endpoint `audited C2-supported reads/M` and for all detection floors are in `tool_audited_signal_auc_comparison.5tools.tsv`.

### C2 Specificity

To test directly whether the strict evidence reflects C2 or Fna signal in general, each tool's `strict-C2-best reads / input reads` is explained in each cohort by both the C2 and the C1 abundance of the original study, adjusting for CRC condition, sex, age, BMI and source sequencing depth. Variables are transformed and standardised within each cohort, and the partial standardised betas are then combined.

| Tool | Published abundance predictor | Three-cohort partial beta | 95% CI | P value |
| --- | --- | ---: | ---: | ---: |
| Chimera | Fna C2 | **0.890** | **0.764 to 1.016** | **1.80e-43** |
| Chimera | Fna C1 | 0.036 | -0.045 to 0.118 | 0.379 |
| Metabuli | Fna C2 | 0.861 | 0.720 to 1.001 | 2.40e-33 |
| Metabuli | Fna C1 | 0.025 | -0.026 to 0.075 | 0.336 |
| ganon | Fna C2 | 0.828 | 0.660 to 0.996 | 4.64e-22 |
| ganon | Fna C1 | 0.024 | -0.038 to 0.086 | 0.445 |

The strict evidence of all three tools mainly follows C2, and none of the C1 betas is significant; Chimera has the highest C2 partial beta.

### CRC Association from Audited Read-level Evidence

To test whether the tool signals recover the CRC association of the original study, this analysis compares CRC with control directly across all 760 samples, without using the paper C2-positive / zero-Fna labels to define the outcome.
In each cohort, the arcsine-square-root transformed endpoint is fitted with the covariate model of the original study (CRC condition, sex, age, BMI and source sequencing depth); the adjusted standardised effects of the three cohorts are then combined with a Paule-Mandel random-effects model.

The read-level endpoint is `audited strict-C2-best reads/M`.
sylph outputs sample-level abundance estimates without per-read assignments and therefore does not apply to this read-level audit.

| Tool / evidence | Three-cohort effect | 95% CI | P value | I2 |
| --- | ---: | ---: | ---: | ---: |
| Published Fna C2 abundance | 0.384 | 0.240 to 0.528 | 1.77e-7 | 0.0% |
| Chimera strict-C2-best reads/M | **0.410** | **0.265 to 0.554** | **2.64e-8** | **0.0%** |
| Metabuli strict-C2-best reads/M | 0.366 | 0.222 to 0.510 | 6.48e-7 | 0.0% |
| Centrifuger strict-C2-best reads/M | 0.214 | 0.071 to 0.357 | 0.00344 | 0.0% |
| ganon strict-C2-best reads/M | 0.201 | 0.057 to 0.344 | 0.00605 | 0.0% |
| Kraken2 LF01 strict-C2-best reads/M | 0.097 | -0.046 to 0.240 | 0.183 | 0.0% |

The three-cohort effect of Chimera's strict evidence is the strongest of the tools tested; its point estimate is above the effect of the published Fna C2 abundance itself, and there is no heterogeneity between the three cohorts.
Chimera's adjusted strict-read effect is positive and significant in each of the three cohorts:

| Cohort | Published Fna C2 effect | Chimera strict-read effect | Chimera P value | Metabuli strict-read effect | Metabuli P value |
| --- | ---: | ---: | ---: | ---: | ---: |
| Yu | 0.437 | **0.581** | 0.00174 | 0.460 | 0.0124 |
| Wirbel | 0.438 | **0.504** | 0.00612 | 0.481 | 0.00870 |
| Yachida | 0.358 | **0.346** | 1.13e-4 | 0.315 | 4.37e-4 |

### CRC Association from Native Signals

The same model applied to the native sample-level C2 signal of each tool:

| Tool / evidence | Signal unit | Three-cohort effect | 95% CI | P value | I2 |
| --- | --- | ---: | ---: | ---: | ---: |
| Published Fna C2 abundance | percent | 0.384 | 0.240 to 0.528 | 1.77e-7 | 0.0% |
| Chimera | reads/M | **0.314** | **0.170 to 0.458** | **1.85e-5** | 0.1% |
| sylph | sequence abundance % | 0.288 | 0.144 to 0.432 | 8.40e-5 | 0.0% |
| Centrifuger | reads/M | 0.207 | 0.063 to 0.350 | 0.00467 | 0.0% |
| Kraken2 LF01 | reads/M | 0.162 | 0.019 to 0.305 | 0.0263 | 0.0% |
| Metabuli | reads/M | 0.133 | -0.013 to 0.279 | 0.0744 | 1.6% |
| ganon | reads/M | 0.131 | -0.012 to 0.274 | 0.0721 | 0.0% |

At the level of native outputs, Chimera has the strongest CRC association of the tools tested. The native signals of Metabuli and ganon do not reach significance; after the orthogonal audit, Metabuli's strict-read effect rises to `0.366`, indicating that its main loss comes from non-specific background in its raw C2 calls, and ganon's is `0.201`.

### Reference and Direct-alignment Robustness

The reference sensitivity analysis uses a pre-specified 60-sample subset whose selection reads no tool result: in each cohort, 10 paper C2-positive samples are chosen at even spacing over their abundance range and matched with 10 zero-Fna samples of the same disease/control composition, giving 30 positive and 30 zero-Fna samples.

The audit reference removes exact duplicate sequences only within the same clade, and keeps identical sequences across clades to represent genuine ambiguity.
Of the original 1,196 reference records, 63 within-clade exact duplicates are removed and 1,133 are kept.
Each read obtains its best alignment within each clade separately, and the `(nmatch, alignment length)` scores are then compared across clades; this design keeps the secondary-hit cap of a combined reference from deciding which clades are visible.
Alignments use minimap2 2.24-r1122 `-x sr` with 128 threads.

The primary robustness endpoint is the all-positive `strict-C2-best reads/M AUC`:

| Method | AUC | Bootstrap 95% CI |
| --- | ---: | ---: |
| Chimera candidate audit, original combined panel | 0.892 | 0.789 to 0.973 |
| Chimera candidate audit, exact-dedup per-clade panel | **0.889** | **0.788 to 0.971** |
| Centrifuger candidate audit, exact-dedup per-clade panel | 0.823 | 0.710 to 0.914 |
| Kraken2 LF01 candidate audit, exact-dedup per-clade panel | 0.745 | 0.613 to 0.858 |
| Direct minimap2 alignment of all reads, exact-dedup per-clade panel | 0.792 | 0.671 to 0.896 |

At `floor >=500`, Chimera ranks all 12 positive samples above the 30 zero-Fna samples under both audit references (AUC `1.000`); with the exact-dedup per-clade panel, Centrifuger reaches `0.964`, direct all-read alignment `0.949` and Kraken2 LF01 `0.918`.
Exact deduplication and per-clade scoring leave Chimera's result essentially unchanged (`0.892` to `0.889`).
Direct all-read alignment shows that high-signal samples can themselves be recovered by targeted alignment, but on the 60-sample subset covering the published abundance range, Chimera's candidate evidence gives a clearer positive/zero-Fna separation.

### Result Files

Paths below are relative to `results/real/`.

| Path | Content |
| --- | --- |
| `fna_c2_crc3_head3m/sample_manifest.tsv` | 760-sample cohort, phenotype, accession, paper abundance, covariate and input-depth manifest |
| `fna_c2_crc3_head3m/sample_level_signals.tsv` | 760 samples x six tools, long-format native Fna C2 / C1 / non-C1/C2 / other signals for plotting |
| `fna_c2_crc3_head3m/read_audit_sample_metrics.tsv` | 760 samples x five read classifiers, candidate fate counts, rates and reads/M |
| `fna_c2_crc3_head3m/reference_manifest.tsv` | frozen species-panel reference: accession, NCBI taxid, species, signal category, source, sequence count, bases, SHA256 and species-label evidence |
| `fna_c2_crc3_head3m/species_panel_dropped.tsv` | genomes removed from the species panel as duplicate strains, with reason |
| `fna_c2_crc3_head3m/audit_reference_manifest.tsv` | frozen paper-clade audit reference: accession, clade, sequence count, bases and SHA256 |
| `fna_c2_crc3_head3m/sample_level_overall_metrics.tsv` | three-cohort overall sample-level C2 signal metrics |
| `fna_c2_crc3_head3m/sample_level_detection_floor_auc.tsv` | three-cohort sample-level detection-floor AUC |
| `fna_c2_crc3_head3m/tool_read_identity_audit/tool_audited_signal_auc_comparison.5tools.tsv` | five-tool read-level audit AUC, bootstrap CI and paired AUC differences |
| `fna_c2_crc3_head3m/tool_read_identity_audit/tool_audited_signal_auc_comparison.4tools.tsv` | the same comparison for Chimera, Centrifuger, Kraken2 LF01 and Metabuli |
| `fna_c2_crc3_head3m/tool_read_identity_audit/tool_audited_signal_auc_comparison.3tools.tsv` | the same comparison for Chimera, Centrifuger and Kraken2 LF01 |
| `fna_c2_crc3_head3m/c2_specificity_partial_association.tsv` | Chimera strict C2 evidence partial association with published C2 and C1 abundance |
| `fna_c2_crc3_head3m/metabuli_c2_specificity_partial_association.tsv` | Metabuli strict C2 evidence partial association with published C2 and C1 abundance |
| `fna_c2_crc3_head3m/ganon_c2_specificity_partial_association.tsv` | ganon strict C2 evidence partial association with published C2 and C1 abundance |
| `fna_c2_crc3_head3m/crc_association_cohort.tsv` | cohort-level CRC association of audited evidence (strict-C2-best reads/M and call precision), bootstrap AUC and adjusted effects |
| `fna_c2_crc3_head3m/crc_association_meta.tsv` | three-cohort random-effects meta-analysis of audited evidence |
| `fna_c2_crc3_head3m/crc_sample_signal_association_cohort.tsv` | cohort-level CRC association for all six tools' native C2 signals |
| `fna_c2_crc3_head3m/crc_sample_signal_association_meta.tsv` | three-cohort native-signal meta-analysis |
| `fna_c2_crc3_head3m/robustness_sample_manifest.tsv` | pre-specified 60-sample reference-sensitivity subset |
| `fna_c2_crc3_head3m/robustness_sensitivity_auc.tsv` | original audit, exact-dedup per-clade audit and all-read direct-alignment AUC |
| `fna_c2_crc3_head3m/reference_panel_exact_dedup_summary.tsv` | per-clade exact-dedup accounting |

### Summary

Across three CRC cohorts and a species-level *Fusobacterium* reference panel, Chimera gives the clearest separation of paper C2-positive from zero-Fna samples on its native output (AUC `0.865`; next best `0.761`) and the lowest background among read classifiers.
On the primary audited endpoint, Chimera reaches `0.897` for all positives and is significantly higher than every comparator at the `>=500` and `>=1000` detection floors.
Chimera's strict evidence follows published C2 rather than C1 abundance and gives the strongest CRC association of all tools (`0.410`; published abundance `0.384`), with no between-cohort heterogeneity.
Exact reference deduplication and per-clade scoring leave Chimera's separation essentially unchanged, and it exceeds direct all-read alignment on the pre-specified 60-sample subset.
