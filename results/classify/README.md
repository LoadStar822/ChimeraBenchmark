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
| centrifuger | cami_refseq | 100 | 142.94611 | 0.908669 | 0.494251 | 0.828105 | 0.414516 | 0.552135 |
| chimera | cami_refseq | 100 | 79.28957 | 0.908669 | 0.900056 | 0.947235 | 0.858027 | 0.900361 |
| ganon | cami_refseq | 100 | 77.28125 | 0.908669 | 0.000196 | 0.642026 | 0.000133 | 0.000266 |
| kraken2 | cami_refseq | 100 | 120.933678 | 0.908669 | 0.400717 | 0.66243 | 0.268262 | 0.381678 |
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
| centrifuger | cami_refseq | 100 | 141.386414 | 0.921764 | 0.852284 | 0.958859 | 0.824495 | 0.886068 |
| chimera | cami_refseq | 100 | 77.162563 | 0.921764 | 0.970719 | 0.95449 | 0.930048 | 0.942092 |
| ganon | cami_refseq | 100 | 76.347656 | 0.921764 | 0.993076 | 0.85283 | 0.849593 | 0.851206 |
| kraken2 | cami_refseq | 100 | 121.849064 | 0.921764 | 0.614647 | 0.85232 | 0.529333 | 0.652511 |
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
| centrifuger | cami_refseq | 10 | 142.233231 | 0.825607 | 0.932399 | 0.911895 | 0.861844 | 0.886156 |
| chimera | cami_refseq | 10 | 78.540939 | 0.825607 | 0.965955 | 0.91668 | 0.894691 | 0.90555 |
| ganon | cami_refseq | 10 | 77.412407 | 0.825607 | 0.934137 | 0.920582 | 0.888068 | 0.904029 |
| kraken2 | cami_refseq | 10 | 122.501129 | 0.825607 | 0.939534 | 0.908604 | 0.860293 | 0.883782 |
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
| centrifuger | cami_refseq | 10 | 141.368347 | 0.829122 | 0.839214 | 0.93136 | 0.819107 | 0.871624 |
| chimera | cami_refseq | 10 | 77.784492 | 0.829122 | 0.956253 | 0.921654 | 0.902307 | 0.911876 |
| ganon | cami_refseq | 10 | 76.571072 | 0.829122 | 0.946583 | 0.919617 | 0.897299 | 0.908319 |
| kraken2 | cami_refseq | 10 | 122.386452 | 0.829122 | 0.833592 | 0.927812 | 0.806868 | 0.863114 |
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
| centrifuger | cami_refseq | 19 | 141.390343 | 1 | 0.598278 | 0.69902 | 0.4365 | 0.526972 |
| chimera | cami_refseq | 19 | 85.92313 | 1 | 0.796233 | 0.818395 | 0.65905 | 0.724692 |
| ganon | cami_refseq | 19 | 77.577396 | 1 | 0.900864 | 0.714675 | 0.64728 | 0.676073 |
| kraken2 | cami_refseq | 19 | 123.761414 | 1 | 0.576879 | 0.696395 | 0.420795 | 0.513281 |
| metabuli | cami_refseq | 19 | 131.255684 | 1 | 0.597423 | 0.728294 | 0.447543 | 0.543185 |

