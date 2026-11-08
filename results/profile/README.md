# Profile Results

For Centrifuger, the tables below report native EM abundance converted to sequence abundance; the read-count comparison is in [`centrifuger_native_comparison.tsv`](centrifuger_native_comparison.tsv) and [`centrifuger_native_sample_metrics.tsv`](centrifuger_native_sample_metrics.tsv).

## Dataset: atcc-hifi

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 120.819668 | 1 | 0.092166 | 0.247183 | 0.283279 |
| centrifuger | cami_refseq | 143.40004 | 1 | 0.208333 | 0.44565 | 0.456659 |
| chimera | cami_refseq | 75.22604 | 1 | 0.952381 | 0.237973 | 0.310359 |
| ganon | cami_refseq | 79.523438 | 1 | 0.487805 | 0.281363 | 0.286109 |
| metabuli | cami_refseq | 128.261719 | 1 | 0.186916 | 0.828272 | 0.723179 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 90.06823 | 0.9 | 0.9 | 0.390159 | 0.372715 |
| sylph | cami_refseq | 22.291855 | 1 | 0.952381 | 0.246766 | 0.43039 |

## Dataset: atcc-illumina

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 121.08617 | 1 | 0.065574 | 0.231177 | 0.115197 |
| centrifuger | cami_refseq | 141.34565 | 1 | 0.031696 | 0.239537 | 0.186683 |
| chimera | cami_refseq | 70.555012 | 1 | 0.952381 | 0.049511 | 0.05169 |
| ganon | cami_refseq | 76.368294 | 1 | 0.065789 | 0.126784 | 0.06825 |
| metabuli | cami_refseq | 128.401966 | 1 | 0.017331 | 0.751987 | 0.581784 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 28.784367 | 0.8 | 0.842105 | 0.168185 | 0.057414 |
| sylph | cami_refseq | 22.247715 | 1 | 1 | 0.053958 | 0.143199 |

## Dataset: cami-strain-madness-long-per-sample

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum; profile metrics are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 100 | 121.117985 | 0.895 | 0.015363 | 0.419096 | 0.16359 |
| centrifuger | cami_refseq | 100 | 142.891922 | 0.9945 | 0.001671 | 0.371746 | 0.170612 |
| chimera | cami_refseq | 100 | 79.28957 | 0.305 | 0.505946 | 0.169846 | 0.07124 |
| ganon | cami_refseq | 100 | 77.28125 | 0.0035 | 0.07 | 1.120811 | 2.672472 |
| metabuli | cami_refseq | 100 | 118.449013 | 0.963 | 0.004169 | 0.327695 | 0.149245 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 100 | 90.175266 | 0 | 0 | 1 | 2.834096 |
| sylph | cami_refseq | 100 | 22.258057 | 0.2275 | 0.948381 | 0.283585 | 0.207443 |

## Dataset: cami-strain-madness-short-per-sample

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum; profile metrics are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 100 | 121.926395 | 0.9565 | 0.016331 | 0.331526 | 0.07053 |
| centrifuger | cami_refseq | 100 | 141.38662 | 0.924421 | 0.144443 | 0.55759 | 0.352078 |
| chimera | cami_refseq | 100 | 77.162563 | 0.569789 | 0.622851 | 0.146144 | 0.047783 |
| ganon | cami_refseq | 100 | 76.347656 | 0.926974 | 0.052344 | 0.363653 | 0.178743 |
| metabuli | cami_refseq | 100 | 128.33522 | 0.999 | 0.012111 | 0.165147 | 0.118899 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 100 | 28.922417 | 0.543789 | 0.62804 | 0.173151 | 0.162584 |
| sylph | cami_refseq | 100 | 22.390266 | 0.547763 | 0.693241 | 0.198332 | 0.128067 |

## Dataset: cami2-marine-long-per-sample

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum; profile metrics are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 10 | 122.735481 | 0.864753 | 0.542685 | 1.34767 | 0.803978 |
| centrifuger | cami_refseq | 10 | 142.233246 | 0.924863 | 0.663819 | 1.341351 | 0.831644 |
| chimera | cami_refseq | 10 | 78.540939 | 0.917978 | 0.802806 | 0.963483 | 0.604775 |
| ganon | cami_refseq | 10 | 77.412407 | 0.903413 | 0.566154 | 1.336373 | 0.815292 |
| metabuli | cami_refseq | 10 | 128.09882 | 0.935011 | 0.216471 | 1.351363 | 0.79947 |
| sylph | cami_refseq | 10 | 23.449219 | 0.871758 | 0.782011 | 1.454061 | 0.885361 |

## Dataset: cami2-marine-short-per-sample

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum; profile metrics are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 10 | 122.414375 | 0.899133 | 0.254249 | 1.338703 | 0.808218 |
| centrifuger | cami_refseq | 10 | 141.456322 | 0.93863 | 0.412803 | 1.325398 | 0.820806 |
| chimera | cami_refseq | 10 | 77.784492 | 0.908958 | 0.773893 | 0.982463 | 0.616264 |
| ganon | cami_refseq | 10 | 76.571072 | 0.920383 | 0.358269 | 1.330662 | 0.817658 |
| metabuli | cami_refseq | 10 | 128.161247 | 0.947726 | 0.061379 | 1.342244 | 0.807309 |
| sylph | cami_refseq | 10 | 23.933979 | 0.857123 | 0.785258 | 1.454177 | 0.881671 |

## Dataset: prjna637878-supported19-single-read

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.
In collection aggregate rows, `Samples` is the number of completed samples; `Max RSS` is the maximum; profile metrics are arithmetic means over samples.

| Tool | DB | Samples | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 19 | 123.798409 | 0.991496 | 0.002046 | 1.595956 | 1.149855 |
| centrifuger | cami_refseq | 19 | 141.35368 | 0.995544 | 0.010039 | 1.422155 | 1.016409 |
| chimera | cami_refseq | 19 | 85.92313 | 0.744754 | 0.15387 | 1.391404 | 1.040058 |
| ganon | cami_refseq | 19 | 77.577396 | 0.989901 | 0.008714 | 1.37134 | 1.023276 |
| metabuli | cami_refseq | 19 | 131.255684 | 0.998051 | 0.001019 | 1.579523 | 1.201258 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 19 | 28.91674 | 0.803589 | 0.124637 | 1.369741 | 1.008034 |
| sylph | cami_refseq | 19 | 24.462078 | 0.73473 | 0.101058 | 1.42739 | 1.076882 |

## Dataset: zymo-gridion-even

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 121.726906 | 0.8 | 0.003945 | 0.401974 | 0.586534 |
| centrifuger | cami_refseq | 142.315861 | 0.8 | 0.000573 | 1.476191 | 1.545048 |
| chimera | cami_refseq | 83.719982 | 0.8 | 1 | 0.198842 | 0.305925 |
| ganon | cami_refseq | 77.804688 | 0.5 | 0.009823 | 1.774869 | 1.668889 |
| metabuli | cami_refseq | 128.17828 | 0.8 | 0.0013 | 0.361536 | 0.538383 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 90.422035 | 0.9 | 1 | 0.221086 | 0.320304 |
| sylph | cami_refseq | 23.300758 | 0.8 | 1 | 0.283743 | 0.437233 |

## Dataset: zymo-gridion-log

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 121.515228 | 0.8 | 0.005768 | 0.091252 | 0.129307 |
| centrifuger | cami_refseq | 142.419666 | 0.8 | 0.000571 | 1.603808 | 2.023259 |
| chimera | cami_refseq | 83.07423 | 0.5 | 1 | 0.025908 | 0.046927 |
| ganon | cami_refseq | 77.966797 | 0.4 | 0.008316 | 1.987853 | 2.189865 |
| metabuli | cami_refseq | 101.861244 | 0.8 | 0.001805 | 0.09416 | 0.138477 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 90.045082 | 0.6 | 1 | 0.018468 | 0.025851 |
| sylph | cami_refseq | 23.306244 | 0.6 | 1 | 0.059437 | 0.100599 |

## Dataset: zymo-promethion-even

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 124.336891 | 0.8 | 0.000523 | 0.398318 | 0.5681 |
| centrifuger | cami_refseq | 142.282619 | 0.8 | 0.002884 | 1.897969 | 1.489213 |
| chimera | cami_refseq | 86.130608 | 0.8 | 1 | 0.187029 | 0.279039 |
| ganon | cami_refseq | 77.810711 | 0.8 | 0.001089 | 1.744334 | 1.613646 |
| metabuli | cami_refseq | 128.365852 | 0.8 | 0.000477 | 0.347398 | 0.517651 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 90.667282 | 0.9 | 1 | 0.207385 | 0.29712 |
| sylph | cami_refseq | 34.440399 | 0.8 | 0.888889 | 0.261798 | 0.405766 |

## Dataset: zymo-promethion-log

### Abundance Metrics

Profiles are scored with the OPAL core metrics; Completeness, Purity and L1 Norm are reported at the `species` rank.
Profiles are compared as sequence abundance.
Only results of the current profile evaluation version are included; results of earlier versions are not in the public tables.
The cross-dataset summary is `profile/summary.tsv`; this page breaks down the results of each tool by dataset.
Weighted UniFrac is a global distance over the whole taxonomy tree; higher Completeness and Purity are better, lower L1 Norm and Weighted UniFrac are better.

| Tool | DB | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Weighted UniFrac |
| --- | --- | --- | --- | --- | --- | --- |
| bracken | cami_refseq | 124.677814 | 0.8 | 0.000559 | 0.090437 | 0.11802 |
| centrifuger | cami_refseq | 142.364391 | 0.8 | 0.003103 | 1.922887 | 1.958933 |
| chimera | cami_refseq | 85.186356 | 0.6 | 1 | 0.025001 | 0.045722 |
| ganon | cami_refseq | 78.451641 | 0.7 | 0.001121 | 1.992073 | 2.240169 |
| metabuli | cami_refseq | 104.923595 | 0.8 | 0.000511 | 0.092679 | 0.135412 |
| metaphlan | mpa_vJan25_CHOCOPhlAnSGB_202503 | 89.940083 | 0.8 | 1 | 0.014503 | 0.011036 |
| sylph | cami_refseq | 34.480919 | 0.7 | 1 | 0.050735 | 0.089043 |
