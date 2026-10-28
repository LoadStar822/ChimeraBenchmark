# ChimeraBenchmark

ChimeraBenchmark is the metagenomic taxonomic classification benchmark used to evaluate Chimera in the accompanying paper. It contains the dataset descriptions, reference collection sizes, evaluation tasks, result tables and the code that produces them.

If you come here from the paper, read `results/README.md` first, then the build, classify and profile result tables.

## Contents

- [Benchmark design](#benchmark-design)
- [Data and results](#data-and-results)
- [Reproduction](#reproduction)
- [Repository layout](#repository-layout)
- [Scope](#scope)
- [License](#license)

## Benchmark design

The benchmark covers three tasks:

- Database build: the cost of building each tool's database from the same reference collection.
- Per-read classification: the accuracy of the assignment of each read or contig at the target rank.
- Abundance profiling: the estimated species composition of each sample.

Results are reported at the `species` rank. Per-read results are scored against per-read truth mappings, and profiles against the composition truth of each sample.

## Data and results

- `results/README.md`: datasets, reference collections, evaluation rules and metrics.
- `results/builds/summary.tsv`, `results/classify/summary.tsv`, `results/classify/sample_metrics.tsv`, `results/profile/summary.tsv`: the result tables used for statistics and figures.
- `results/paper_run_manifest.tsv`: the runs behind the result tables, with tool versions, parameter policies and completion status.
- `results/real/README.md`: design and results of the real-cohort experiment (the *Fusobacterium animalis* signal in colorectal cancer cohorts); its tables are in `results/real/fna_c2_crc3_head3m/`.
- `results/builds/README.md`: database build cost and the RefSeq scaling experiment. `results/classify/README.md` and `results/profile/README.md`: the result tables broken down by dataset.
- Resources that are too large for this repository, such as the NCBI taxonomy snapshot and the derived truth tables, are archived at Zenodo ([doi:10.5281/zenodo.23124825](https://doi.org/10.5281/zenodo.23124825)) and unpack under `resources/`.
- The Chimera database used for all results can be downloaded pre-built; see [Pre-built Database](https://github.com/LoadStar822/Chimera#pre-built-database) in the Chimera README.

Results come in three layers: a run is one execution in an experiment directory; an evaluation is the set of metrics scored for that run; the result tables include only runs that were explicitly selected and passed validation, and these runs are listed in `results/paper_run_manifest.tsv`. Running and evaluating never modify the result tables, and generating the result tables never reruns a tool.

## Reproduction

List the available commands:

```bash
python -m chimera_bench.cli -h
```

Regenerate the data resource catalog (writes only to `resources/reports/`):

```bash
python -m chimera_bench.cli catalog --config configs --resources-root resources
```

The configurations in `configs/` record the input, truth and database paths used on the benchmark server; point them to local copies before running.

## Repository layout

- `configs/`: dataset, reference database and tool run configurations.
- `chimera_bench/`: code for running, summarising and evaluating.
- `results/`: result tables, per-dataset breakdowns and experiment notes.
- `resources/` (not tracked; see the Zenodo archive above): the taxonomy snapshot, input lists, derived truth tables, the reference genomes of the real-cohort experiment and the data resource catalog.

## Scope

- Apart from thread count and output paths, every tool runs with its default parameters. Settings chosen by input or output type, such as Metabuli's `--seq-mode` and the `--binning` preset used for ganon's per-read results, are listed in `results/README.md`.
- Peak memory of multi-step tools is measured over the whole pipeline, for example `kraken2 -> bracken`.
- All evaluations use the NCBI taxonomy snapshot `resources/taxonomy/ncbi_20260408/`.
- Per-read F1 should be read together with the `Truth Mapped Rate`: when the mapped rate is low, the result describes only the mappable subset.

## License

The code in this repository is released under the MIT License (see `LICENSE`).
