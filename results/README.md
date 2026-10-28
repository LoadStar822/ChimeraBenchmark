# Benchmark Data and Results

This document describes the data, reference collections, evaluation tasks, metrics and result tables of ChimeraBenchmark.
It is written for readers who want to check data sources, evaluation rules and results.

## Contents

- [Data and result files](#data-and-result-files)
- [Using the results](#using-the-results)
- [Commands](#commands)
- [Reference Database Summary](#reference-database-summary)
- [Benchmark Dataset Summary](#benchmark-dataset-summary)
- [Evaluation tasks and rules](#evaluation-tasks-and-rules)
- [Dataset-specific rules](#dataset-specific-rules)
- [Result summary](#result-summary)
- [Real-cohort experiment](#real-cohort-experiment)
- [Metrics](#metrics)
- [Tools](#tools)
- [Software versions](#software-versions)

## Data and result files

- `builds/summary.tsv`, `classify/summary.tsv`, `classify/sample_metrics.tsv`, `profile/summary.tsv`: the result tables. `classify/README.md` and `profile/README.md` break them down by dataset.
- `paper_run_manifest.tsv`: the selected runs, with tool versions, parameter policies, completion status and result files.
- `results/real/README.md`: design and interpretation of the real-cohort experiment; its tables are in `results/real/fna_c2_crc3_head3m/`.
- `results/builds/README.md`: database build cost and the RefSeq scaling experiment.
- `resources/reports/db_catalog.tsv`: TSV version of the Reference Database Summary.
- `resources/reports/dataset_catalog.tsv`: TSV version of the Benchmark Dataset Summary.

`resources/` is not tracked in this repository; its contents are archived at Zenodo ([doi:10.5281/zenodo.23124825](https://doi.org/10.5281/zenodo.23124825)).

Results are reported at the `species` rank.

Runs, evaluations and result tables are three different objects: a run belongs to one explicit experiment directory; an evaluation records the scores of that run; the result tables include only runs that were explicitly selected and passed validation. `paper_run_manifest.tsv` lists these runs with their tool versions, parameter policies and purpose (main or supplementary experiment).
Ad hoc experiments and metric recomputation never update the result tables. This document and the other notes are maintained by hand; the catalog command generates only `resources/reports/catalog.md` and the accompanying TSV files.

## Using the results

- All tools in the default benchmark run with the same thread count, `threads=32`. The thread settings of the Fna C2 real-cohort experiment are given in `results/real/README.md`.
- Apart from thread count and output paths, every tool runs with its default parameters. Settings that a tool's documentation ties to the input type follow the dataset: Metabuli's `--seq-mode`, MetaPhlAn 4's `--long_reads` and Bracken's read-length distribution. ganon's per-read results use the `--binning` preset that its documentation gives for read binning. See [Tools](#tools); the parameter policy of each entry is recorded in the `parameter_policy` column of `paper_run_manifest.tsv`.
- All classifiers use databases built from the same reference collection (CAMI RefSeq); MetaPhlAn 4 uses its own marker gene database.
- All evaluations use the NCBI taxonomy snapshot `resources/taxonomy/ncbi_20260408/`.
- CAMI II Strain Madness (100 long-read and 100 short-read samples) and CAMI II Marine (10 long-read and 10 short-read samples) are run sample by sample. Accuracy metrics in the summary tables are arithmetic means over samples, and peak memory is the maximum over samples; per-sample results are in `classify/sample_metrics.tsv`.

## Commands

Run the following commands from the ChimeraBenchmark root. `EXP` is a configuration name in `configs/experiments/`. Choose a new `RUN` for every new experiment and never point it at an aggregated result directory.

```bash
EXP=chimera
RUN=results/experiments/chimera-trial
python -m chimera_bench.cli run --exp "$EXP" --config configs --runs "$RUN"
python -m chimera_bench.cli recompute --exp "$EXP" --config configs --runs "$RUN"
python -m chimera_bench.cli report --runs "$RUN" --out "$RUN/summary.tsv"
```

- `run` and `build` require an explicit `--runs`. The root directory is tagged by `experiment.json` with the experiment identity and a digest of the effective parameters. Different experiments or parameters cannot share a root, and an existing run directory cannot be reused.
- Classification results, native profiles and logs are kept in the same experiment directory; a separate `--profile` root is no longer accepted.
- Build outputs go into the directory of that build and never overwrite the shared database that a configuration points to. The `db` of later classification configurations should be the actual database path recorded under `outputs` in the build's `meta.json`.
- `recompute` does not rerun the tool. It replaces the metrics of the selected completed runs and records the evaluation configuration it used. A failed first scoring can be recovered by recomputing; runs that are still executing or being scored are not recomputed.
- `report` exports only to `.tsv` files directly inside the experiment root and cannot overwrite READMEs, identity tags or raw outputs in subdirectories. If no complete result matches, it refuses to export and keeps any existing file.

## Reference Database Summary

This table summarises the size of the reference collections used by the benchmark.
It describes the reference collections themselves, not the database directory built by any tool.
The statistics come from the reference target TSV and NCBI assembly summary metadata.

| Dataset Name | Total Size (GB) | Total Sequences | Base Pairs (bp) | Assemblies | Species Count | Source |
| --- | --- | --- | --- | --- | --- | --- |
| CAMI RefSeq | 173.0 | 16,738,778 | 588,675,466,098 | 141,677 | 21,075 | NCBI RefSeq genomic, CAMI snapshot |
| RefSeq Complete | 55.3 | 119,431 | 189,355,217,616 | 58,666 | 25,126 | NCBI RefSeq complete genomes |

## Benchmark Dataset Summary

This table summarises the input files of the 11 datasets in the result tables.
The statistics come from `configs/datasets/*.yaml`, from scanning the input sequence files, or from metadata shipped with each dataset.

| Dataset Name | Total Size (GB) | Samples | Input Type | Reads / Contigs | Base Pairs (bp) | Mean Length (bp) | N50 (bp) | GC (%) | Q30 (%) | Truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CAMI II Marine (long-read contigs) | 10.8 | 10 | contig FASTA | 1,622,375 | 10,587,029,583 | 6,526 | 9,364 | 48.12 | — | per-contig mapping; profile weighted by reads per contig |
| CAMI II Marine (short-read contigs) | 10.8 | 10 | contig FASTA | 18,057,427 | 10,389,109,998 | 575 | 770 | 48.11 | — | per-contig mapping; profile weighted by reads per contig |
| CAMI II Strain Madness (long reads) | 401.1 | 100 | single FASTQ | 69,647,201 | — | — | — | — | — | per-read mapping; profile by read bases |
| CAMI II Strain Madness (short reads) | 424.2 | 100 | paired FASTQ | 1,330,822,164 | 199,623,324,600 | 150 | 150 | — | — | per-read mapping; profile by read count (2×150 bp) |
| ATCC MSA-1003 (Illumina) | 3.3 | 1 | paired FASTQ | 10,038,314 | 1,254,789,250 | 125 | 125 | 48.91 | 92.54 | profile only (DNA composition) |
| ATCC MSA-1003 (PacBio HiFi) | 41.3 | 1 | single FASTQ | 2,419,037 | 20,543,987,923 | 8,493 | 9,198 | 53.75 | 96.70 | profile only (DNA composition) |
| ZymoBIOMICS Even (GridION) | 28.5 | 1 | single FASTQ | 3,491,078 | 14,007,156,825 | 4,012 | 5,213 | 46.02 | 0.00 | profile only (theoretical DNA composition) |
| ZymoBIOMICS Log (GridION) | 32.6 | 1 | single FASTQ | 3,667,007 | 16,032,264,247 | 4,372 | 5,290 | 40.94 | 0.00 | profile only (theoretical DNA composition) |
| ZymoBIOMICS Even (PromethION) | 298.4 | 1 | single FASTQ | 36,527,376 | 146,291,709,850 | 4,005 | 5,288 | 46.38 | 0.00 | profile only (theoretical DNA composition) |
| ZymoBIOMICS Log (PromethION) | 301.6 | 1 | single FASTQ | 35,118,078 | 148,028,914,788 | 4,215 | 5,219 | 41.41 | 0.00 | profile only (theoretical DNA composition) |
| PRJNA637878 Supported Fecal Genomes (single-read) | 110.5 | 19 | single FASTQ | 1,321,826,777 | 169,812,107,628 | 128 | — | — | — | mate-level per-read; profile from the same labels |

## Evaluation tasks and rules

Per-read classification evaluates the assignment of each read or contig.
Abundance profiling evaluates the estimated species composition of the whole sample.

For species-level per-read results, both the truth taxonomy ID (taxid) and the predicted taxid are first lifted to the species rank and then compared.
This treats truth labels given at strain, subspecies or isolate level fairly.

Truth and predicted profiles are both compared as sequence abundance (the share of the sample's sequence attributed to each taxon), not as cell counts or genome copies:

- CAMI II Strain Madness truth is summed from the per-read mapping: long reads by the bases of each species' reads (`resources/derived/strain_madness_long_bases/`); short reads all have the same length (2×150 bp) and are counted by reads.
- CAMI II Marine is weighted by the number of source reads of each contig (the `number_reads` column of the mapping), as in the official CAMI rules.
- ATCC MSA-1003 and ZymoBIOMICS use DNA-based compositions; see [Dataset-specific rules](#dataset-specific-rules).
- The PRJNA637878 profile truth counts reads of the same per-read labels.

Predicted abundances are taken from each tool's native output without unit conversion: Chimera reports abundance by bases; Sylph uses the `Sequence_abundance` column; MetaPhlAn 4 uses the estimated read counts (bases for long reads) of `-t rel_ab_w_read_stats`; ganon uses `ganon report --report-type reads`; Bracken, Metabuli and Centrifuger use their default abundance outputs.

## Dataset-specific rules

### CAMI II Marine

The inputs are the gold standard assembly (GSA) contigs provided by CAMI. Per-read metrics are computed per contig, and each contig counts once.
MetaPhlAn 4 targets sequencing reads and is not run on these two contig datasets.

### ATCC MSA-1003

ATCC MSA-1003 is a staggered mock community of 20 species without per-read truth; only profiles are evaluated.
The truth `configs/truth/atcc_mock_dna.tsv` multiplies the genome copy proportions stated for the product by the genome size of each strain and normalises them, giving a DNA-based composition. How the genome sizes were chosen is documented in the comments of that file.

### ZymoBIOMICS

ZymoBIOMICS comprises an even standard (8 bacteria, 2 yeasts) and a log-distributed standard, each sequenced on GridION and PromethION, giving 4 datasets; only profiles are evaluated.
The truth is the theoretical composition by genomic DNA given by the manufacturer for the even standard D6300 (`configs/truth/zymo_d6300.tsv`) and the log-distributed standard D6310 (`configs/truth/zymo_d6310.tsv`).

The two yeasts (*Saccharomyces cerevisiae*, *Cryptococcus neoformans*) have no genomes in the CAMI RefSeq reference collection.
Tools using this collection therefore reach at most `0.8` completeness on ZymoBIOMICS. In the log-distributed standard *S. cerevisiae* makes up 0.89%, which puts a floor of about `0.018` under the species L1 norm of these tools.
MetaPhlAn 4 uses its own marker gene database, which includes *S. cerevisiae*.

### PRJNA637878

PRJNA637878 pairs stool metagenomes with genomes of bacterial isolates cultured from the same samples. The benchmark uses the 19 human stool samples (10 healthy, 5 recurrent *C. difficile* infection, 4 inflammatory bowel disease) for which isolates from the same sample were sequenced.

The per-read truth covers only reads that can be traced to a cultured isolate. The paired reads of each sample were aligned with minimap2 (`-ax sr -N 20`) to all 956 isolate assemblies of the project. A mate is labelled with the species of its best alignment when that alignment spans at least 140 bp with at most 10 bp soft-clipped, all alignments within 5 alignment-score points of the best agree on one species, that species was isolated from the same sample, and the highest mapping quality among these alignments is at least 30. Species labels are the genus and species given for each isolate in Supplementary Table 1 of the original publication. Mates without a label are not scored.

Some sequencing runs (SRA runs) of PRJNA637878 contain unpaired (singleton) reads, visible as `spots_with_mates` smaller than `spots`.
`prjna637878-supported19-single-read` therefore evaluates R1 and R2 as independent single reads (single-read evaluation), and adds a mate tag `__mate1` or `__mate2` to each read identifier when generating the FASTQ and truth tables.
This does not change the paired-end origin of the original experiment; it only fixes that all tools are evaluated on the same single-read input without mate-pair information.

#### Truth correction

A mate's truth label is the species name of the best-matching assembly, which may come from another sample of the project. Errors in those names therefore enter the truth. We checked all 956 assemblies against the CAMI RefSeq reference collection with genome-level average nucleotide identity (ANI, sylph naive ANI):

- A species is accepted for an assembly when more than half of its reference genomes have ANI ≥ 95% to the assembly.
- A label is kept when its species is accepted, when it is the species of the best match, or when it reaches at least 95% ANI within 0.5 points of the best match.
- Otherwise the label is replaced by the accepted species with the highest ANI, preferring formal binomial names. Labels with no accepted species (novel species) are kept.

The audit depends only on the assemblies and the reference collection, not on any tool's output. Among the assemblies whose label changes the species taxid, nine carry a species that is present in the reference collection; their labels are corrected below. Name changes that keep the species taxid (for example *Bacteroides massiliensis* to *Phocaeicola massiliensis*) need no correction, and labels of species absent from the reference collection are left unchanged.

| Assembly | Published species | Corrected species | ANI to published species (%) | ANI to corrected species (%) | Relabelled mates | NCBI record (checked 2026-09-23) |
| --- | --- | --- | ---: | ---: | ---: | --- |
| GCF_015667855.1 | *Blautia wexlerae* | *Fusicatenibacter saccharivorans* | 87.43 | 98.18 | 317,528 | suppressed (unverified source organism) |
| GCF_015668735.1 | *Bacteroides thetaiotaomicron* | *Bacteroides faecis* | 93.55 | 99.25 | 298,355 | *Bacteroides faecis* |
| GCF_015669505.1 | *Lactobacillus rhamnosus* | *Pediococcus acidilactici* | n.d. | 99.78 | 31,961 | *Pediococcus acidilactici* |
| GCF_015668165.1 | *Blautia massiliensis* | *Blautia wexlerae* | 88.06 | 97.80 | 5,974 | renamed *Blautia* sp. |
| GCF_015553015.1 | *Streptococcus oralis* | *Streptococcus halitosis* | 98.76 | 96.19 | 348 | *Streptococcus oralis* |
| GCF_015668095.1 | *Parabacteroides johnsonii* | *Parabacteroides merdae* | 92.99 | 99.89 | 32 | *Parabacteroides merdae* |
| GCF_015667505.1 | *Clostridium tepidum* | *Faecalibacillus intestinalis* | n.d. | 98.78 | 0 | *Faecalibacillus intestinalis* |
| GCF_015670125.1 | *Clostridium colicanis* | *Clostridium thermobutyricum* | n.d. | 98.62 | 0 | *Clostridium thermobutyricum* |
| GCF_015668845.1 | *Paeniclostridium sordellii* | *Paraclostridium bifermentans* | 87.70 | 97.67 | 0 | *Paraclostridium bifermentans* |

ANI values are the best match among the reference genomes of each species; n.d., no match detected. For *S. oralis*, 28 of the 107 reference genomes reach 95% ANI, so the species is not accepted under the majority rule, while the single *S. halitosis* genome is.

Mates carrying one of the nine published names were relabelled as follows. Each mate and its pair were aligned with minimap2 (`-x sr`, secondary alignments kept) to the 956 assemblies, and only assemblies with the same published name as the mate's label were considered. An assembly scores the best alignment score of the mate (matches minus mismatches) plus that of its pair on the same assembly. A mate was relabelled only when all top-scoring assemblies share the same corrected species and none of them keeps its label.

In total 654,198 of 42,090,303 labelled mates (1.55%) in 14 of the 19 samples were relabelled. All tools are scored against the corrected truth; the original and corrected truth tables are both in the Zenodo archive.

## Result summary

Cross-dataset summaries use only explicitly selected complete runs (see `paper_run_manifest.tsv`) and never scan directories for the "latest run". All TSV files and per-dataset tables come from the same set of runs.

- Per-read classification covers CAMI II Strain Madness long/short, CAMI II Marine long/short and `prjna637878-supported19-single-read`, and compares Chimera, Centrifuger, Kraken2, Metabuli and ganon. ATCC and ZymoBIOMICS have no per-read truth.
- Abundance profiling covers all 11 datasets and compares the native abundance outputs of Chimera, Centrifuger, Metabuli, ganon, Bracken, Sylph and MetaPhlAn 4. Abundance estimates for Kraken2 are given by Bracken.

The summary tables are `builds/summary.tsv`, `classify/summary.tsv`, `classify/sample_metrics.tsv` and `profile/summary.tsv`.

## Real-cohort experiment

The Fna C2 real-cohort experiment is in `results/real/`.
On 760 stool metagenomes from three colorectal cancer (CRC) cohorts, it uses a species-level reference panel of 188 reference genomes from 8 *Fusobacterium* species to evaluate how well each tool detects reads of the *Fusobacterium animalis* C2 clade (Fna C2), and how this signal associates with CRC.

`results/real/README.md` summarises the design, metrics, results and accompanying tables.

The detection floor converts the C2 abundance reported in the original paper and the input depth into expected C2 reads, and keeps only C2-positive samples that reach a given floor in the positive group.
This separates "is there a C2 signal at all" from "once the C2 signal is strong enough to detect, does the tool separate the samples".

## Metrics

### Per-read metrics

Per-read metrics measure whether each sequence is assigned to the correct taxon at the target rank.

- `Precision`: of the sequences a tool assigns at the target rank, the fraction that is correct.
- `Recall`: of the sequences that can be evaluated in the truth, the fraction recovered correctly.
- `F1`: the harmonic mean of precision and recall.
- `Truth Mapped Rate`: the fraction of truth sequences that map to the target rank.
- `Pred Mapped Rate`: the fraction of predicted sequences that map to the target rank.

When the mapped rate is clearly below `1`, F1 describes performance on the mappable subset and should not be read alone as performance on the whole dataset.

### Profile metrics

The profile tables use the OPAL core metrics:
- `Completeness`: the fraction of taxa present in the truth that are predicted.
- `Purity`: the fraction of predicted taxa that are present in the truth.
- `L1 Norm`: the L1 distance between the truth and predicted abundance distributions, ranging from `0` to `2`.
- `Weighted UniFrac`: the difference of cumulative abundances along the taxonomy tree.

Higher completeness and purity are better; lower L1 norm and weighted UniFrac are better.

Profile results use only the native abundance outputs of the tools.
A tool without native profile output does not appear in the profile tables.

### Fields

- `Elapsed Seconds`: database build time in seconds.
- `Max RSS (GB)`: peak resident memory, with 1 GB = 1024³ bytes.
- `DB Size`: file size of the database directory (`du --apparent-size`), in GiB.
- `Total Size (GB)`: size of the input or reference files in decimal GB.
- `Q30 (%)`: the fraction of FASTQ bases with quality score at least 30.

## Tools

### Chimera

Per-read results and abundance profiles come from the same `classify` run with default parameters.
`--profile-cami` only requests an additional CAMI-format abundance file and does not change read classification; that file comes from the same estimate as the native abundance output and is by bases.

### Kraken2

Kraken2 provides only per-read classification; its abundance estimates are given by Bracken.

### Bracken

Bracken re-estimates abundance from Kraken2 output and provides no per-read classification.
It therefore appears only in the profile tables.
Its peak memory is measured over the whole `kraken2 -> bracken` pipeline.

The Bracken database provides `database{L}mers.kmer_distrib` for read lengths of 100, 125 and 150 bp.
Following the Bracken manual, datasets with a stated read length use the matching distribution (125 bp for ATCC Illumina; 150 bp for CAMI II Strain Madness short and PRJNA637878), and the other datasets use the default 100 bp.

### Centrifuger

Centrifuger's per-read classification uses the default inference parameters.
Profiles are produced by `centrifuger-quant` from the classification output and exported in CAMI profile format.

### Metabuli

Metabuli uses default parameters. For single-file input, `--seq-mode` follows its documentation: `1` for single-end short reads and `3` for long reads and contigs; paired input uses the default mode.
Profiles are converted from the Metabuli report to CAMI profile format.

### ganon

ganon builds a species-level database on CAMI RefSeq (`ganon build-custom`; species-level targets are the default level of `ganon build`). The build uses `--skip-genome-size`, which only affects abundance reports converted to genome copies and not the sequence abundance used in the evaluation.

Following the ganon documentation, per-read classification and abundance profiling are run separately: per-read classification uses the `--binning` preset optimised for read binning (`--rel-cutoff 0.25 --rel-filter 0 --min-count 0`), and profiling uses `classify` with default parameters followed by `ganon report --report-type reads` for sequence abundance.
Both classifications belong to the same task, and peak memory is measured over the whole task.

On the simulated long reads with high error rates (CAMI II Strain Madness long), ganon classifies almost nothing under its default thresholds, and both its per-read and abundance results are close to zero.

### Sylph

Sylph uses the default sketch parameters and the default `profile`; abundance is taken from the `Sequence_abundance` column.
Sylph provides no per-read classification and appears only in the profile tables.

### MetaPhlAn 4

MetaPhlAn 4 uses its own marker gene database `mpa_vJan25_CHOCOPhlAnSGB_202503` instead of the CAMI RefSeq reference collection, and appears only in the profile tables.
Short reads are aligned with bowtie2 by default; for long-read platforms (PacBio HiFi, Nanopore, CAMI II Strain Madness long) `--long_reads` is added as documented, which switches to minimap2.
The output uses `-t rel_ab_w_read_stats`, and its estimated read counts (bases for long reads) are taken as sequence abundance.

On CAMI II Strain Madness long, MetaPhlAn 4 reports no species.
In long-read mode the default filter requires the sequence divergence of an alignment (the minimap2 `de` tag) to be at most `--max_gcsd 0.10`, while the simulated reads of this dataset diverge more from the marker genes: realigning the first 2,000 reads of sample 0 with MetaPhlAn 4's long-read alignment parameters, 86 reads align to marker genes (204 alignment records) with a median divergence of 0.160 and a 5th percentile of 0.142, and no alignment satisfies `de ≤ 0.10`.

## Software versions

- chimera: 1.8.0
- kraken2: 2.1.3
- bracken: 2.9 (output of `bracken -v`; the bioconda package version is 3.0)
- centrifuger: 1.1.0-r291
- metabuli: 1.2.0
- ganon: 2.1.0
- sylph: 0.8.1
- metaphlan: 4.2.2 (database `mpa_vJan25_CHOCOPhlAnSGB_202503`)
