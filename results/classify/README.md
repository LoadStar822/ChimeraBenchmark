# Classify Results

Auto-generated. Do not edit.

## Dataset: cami-strain-madness-long-per-sample

### Per-read Metrics

Classification is judged at the `species` rank: truth and prediction are both lifted to the species rank and then compared.
For example, when the truth is a species and the prediction is a strain of that species, the `species` column counts it as correct.
Truth labels at `strain/subspecies/isolate` level enter the denominator if they can be lifted to a species; truth labels only at `genus` or higher do not.
Truth Mapped Rate and Pred Mapped Rate are shown alongside; F1 must be read together with the mapped rates and not quoted apart from its denominator.
The cross-dataset summary is `classify/summary.tsv`; this page breaks down the results of each tool by dataset.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum and read counts are summed; rates and accuracies are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 100 | 142.891922 | 0.908669 | 0.493782 | 0.827928 | 0.414038 | 0.551673 |
| chimera | cami_refseq | 100 | 79.28957 | 0.908669 | 0.900056 | 0.947235 | 0.858027 | 0.900361 |
| ganon | cami_refseq | 100 | 77.28125 | 0.908669 | 0.000196 | 0.642026 | 0.000133 | 0.000266 |
| kraken2 | cami_refseq | 100 | 120.923401 | 0.908669 | 0.400184 | 0.663133 | 0.26782 | 0.381345 |
| metabuli | cami_refseq | 100 | 118.449013 | 0.908669 | 0.260807 | 0.807561 | 0.211236 | 0.334714 |

## Dataset: cami-strain-madness-short-per-sample

### Per-read Metrics

Classification is judged at the `species` rank: truth and prediction are both lifted to the species rank and then compared.
For example, when the truth is a species and the prediction is a strain of that species, the `species` column counts it as correct.
Truth labels at `strain/subspecies/isolate` level enter the denominator if they can be lifted to a species; truth labels only at `genus` or higher do not.
Truth Mapped Rate and Pred Mapped Rate are shown alongside; F1 must be read together with the mapped rates and not quoted apart from its denominator.
The cross-dataset summary is `classify/summary.tsv`; this page breaks down the results of each tool by dataset.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum and read counts are summed; rates and accuracies are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 100 | 141.38662 | 0.921764 | 0.85197 | 0.958748 | 0.824423 | 0.885979 |
| chimera | cami_refseq | 100 | 77.162563 | 0.921764 | 0.970719 | 0.95449 | 0.930048 | 0.942092 |
| ganon | cami_refseq | 100 | 76.347656 | 0.921764 | 0.993076 | 0.85283 | 0.849593 | 0.851206 |
| kraken2 | cami_refseq | 100 | 121.870834 | 0.921764 | 0.613312 | 0.853477 | 0.528374 | 0.652116 |
| metabuli | cami_refseq | 100 | 128.33522 | 0.921764 | 0.716911 | 0.926741 | 0.669749 | 0.777093 |

## Dataset: cami2-marine-long-per-sample

### Per-read Metrics

Classification is judged at the `species` rank: truth and prediction are both lifted to the species rank and then compared.
For example, when the truth is a species and the prediction is a strain of that species, the `species` column counts it as correct.
Truth labels at `strain/subspecies/isolate` level enter the denominator if they can be lifted to a species; truth labels only at `genus` or higher do not.
Truth Mapped Rate and Pred Mapped Rate are shown alongside; F1 must be read together with the mapped rates and not quoted apart from its denominator.
The cross-dataset summary is `classify/summary.tsv`; this page breaks down the results of each tool by dataset.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum and read counts are summed; rates and accuracies are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 10 | 142.233246 | 0.825607 | 0.943408 | 0.907173 | 0.865867 | 0.886035 |
| chimera | cami_refseq | 10 | 78.540939 | 0.825607 | 0.965955 | 0.91668 | 0.894691 | 0.90555 |
| ganon | cami_refseq | 10 | 77.412407 | 0.825607 | 0.934137 | 0.920582 | 0.888068 | 0.904029 |
| kraken2 | cami_refseq | 10 | 122.690838 | 0.825607 | 0.94009 | 0.907466 | 0.86154 | 0.883902 |
| metabuli | cami_refseq | 10 | 128.09882 | 0.825607 | 0.898604 | 0.917379 | 0.841806 | 0.877962 |

## Dataset: cami2-marine-short-per-sample

### Per-read Metrics

Classification is judged at the `species` rank: truth and prediction are both lifted to the species rank and then compared.
For example, when the truth is a species and the prediction is a strain of that species, the `species` column counts it as correct.
Truth labels at `strain/subspecies/isolate` level enter the denominator if they can be lifted to a species; truth labels only at `genus` or higher do not.
Truth Mapped Rate and Pred Mapped Rate are shown alongside; F1 must be read together with the mapped rates and not quoted apart from its denominator.
The cross-dataset summary is `classify/summary.tsv`; this page breaks down the results of each tool by dataset.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum and read counts are summed; rates and accuracies are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 10 | 141.456322 | 0.829122 | 0.849435 | 0.92696 | 0.822576 | 0.871645 |
| chimera | cami_refseq | 10 | 77.784492 | 0.829122 | 0.956253 | 0.921654 | 0.902307 | 0.911876 |
| ganon | cami_refseq | 10 | 76.571072 | 0.829122 | 0.946583 | 0.919617 | 0.897299 | 0.908319 |
| kraken2 | cami_refseq | 10 | 122.431698 | 0.829122 | 0.835791 | 0.928602 | 0.80999 | 0.865239 |
| metabuli | cami_refseq | 10 | 128.161247 | 0.829122 | 0.831588 | 0.93014 | 0.804595 | 0.862815 |

## Dataset: prjna637878-supported19-single-read

### Per-read Metrics

Classification is judged at the `species` rank: truth and prediction are both lifted to the species rank and then compared.
For example, when the truth is a species and the prediction is a strain of that species, the `species` column counts it as correct.
Truth labels at `strain/subspecies/isolate` level enter the denominator if they can be lifted to a species; truth labels only at `genus` or higher do not.
Truth Mapped Rate and Pred Mapped Rate are shown alongside; F1 must be read together with the mapped rates and not quoted apart from its denominator.
The cross-dataset summary is `classify/summary.tsv`; this page breaks down the results of each tool by dataset.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum and read counts are summed; rates and accuracies are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Truth Mapped Rate (species) | Pred Mapped Rate (species) | Precision (species) | Recall (species) | F1 (species) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| centrifuger | cami_refseq | 19 | 141.35368 | 1 | 0.600709 | 0.701476 | 0.438901 | 0.529357 |
| chimera | cami_refseq | 19 | 85.92313 | 1 | 0.796233 | 0.818395 | 0.65905 | 0.724692 |
| ganon | cami_refseq | 19 | 77.577396 | 1 | 0.900864 | 0.714675 | 0.64728 | 0.676073 |
| kraken2 | cami_refseq | 19 | 123.732159 | 1 | 0.575178 | 0.694918 | 0.418794 | 0.510976 |
| metabuli | cami_refseq | 19 | 131.255684 | 1 | 0.597423 | 0.728294 | 0.447543 | 0.543185 |

