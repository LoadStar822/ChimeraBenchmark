# Classify Directory README
# 1. Experiment Overview

This directory is dedicated to metagenomic classification experiments, including scripts for running classification tasks, example input data, and result analysis. Similar to the **Build** experiment, this study evaluates the performance of multiple classification tools, including **Chimera** and other alignment-based software.

## Experiment Objective

This experiment aims to assess the performance of different classification tools in metagenomic data classification, focusing on the following aspects:

- **Classification Accuracy**: Evaluating the correctness of classifications across various metagenomic datasets.
- **Processing Efficiency**: Measuring execution time and computational cost at different dataset scales.
- **Resource Consumption**: Analyzing memory usage and storage requirements for classification tasks.

Through this experiment, we seek to compare the performance of **Chimera** with other classification tools, providing data-driven insights for further optimization and tool selection.

## Evaluation Metrics and Formulas

The classification performance is evaluated using the following five metrics:

### 1. Accuracy
Represents the proportion of correctly classified instances. It is calculated as:

$$
Accuracy = \frac{TP + TN}{TP + TN + FP + FN}
$$

### 2. Precision
Measures the correctness of the classifier when predicting positive categories. It is calculated as:

$$
Precision = \frac{TP}{TP + FP}
$$

### 3. Recall (Sensitivity)
Represents the proportion of actual positive instances correctly identified by the classifier. It is calculated as:

$$
Recall = \frac{TP}{TP + FN}
$$

### 4. F1 Score
The harmonic mean of precision and recall, balancing their effects. It is calculated as:

$$
F1 = 2 \times \frac{Precision \times Recall}{Precision + Recall}
$$

### 5. L1 Distance
Measures the deviation between classification results and the true taxonomic distribution. It is calculated as:

$$
L1\_Distance = \sum_{i} | P_i - T_i |
$$

where:
- $P_i$ represents the normalized abundance of the predicted category.
- $T_i$ represents the normalized abundance of the true category.

## Classification Result Definitions

- **TP (True Positive)**: DNA sequences correctly classified into their true species.
- **FN (False Negative)**: DNA sequences that failed to be classified.
- **FP (False Positive)**: DNA sequences misclassified as another species.
- **TN (True Negative)**: DNA sequences correctly classified as non-target species.

These metrics comprehensively evaluate the classification performance of **Chimera** and other tools, considering classification accuracy, misclassification impact, and distributional deviations.

# 2. Experimental Data

## Datasets

This experiment employs four real and simulated datasets from the **CAMI II** project, supplemented by an additional simulated dataset to enhance classification evaluation coverage (see Supplemental Table S3). These datasets encompass diverse sequencing strategies, including **long-read sequencing (average ~3,000 bp)** and **short-read sequencing (2 × 150 bp)**, designed to simulate complex microbial ecosystems such as **marine** and **mouse gut microbiomes**. 

Additionally, the simulated datasets incorporate **technical sequencing errors** and **random insert-size variations** to assess the robustness and adaptability of classification algorithms.

### Dataset Details

| Dataset Name | Total Size (GB) | Samples | Read Length (bp) | Read Length S.D. (bp) | Insert Size Mean (bp) | Insert Size S.D. (bp) |
|-------------|----------------|---------|------------------|----------------------|----------------------|----------------------|
| **CAMI II Marine (long read)** | 10.1 | 10 | Average 3,000 | 1,000 | — | — |
| **CAMI II Marine (short read)** | 11.5 | 10 | 2 × 150 | — | 270 | 20 |
| **CAMI II Toy Mouse Gut (long read)** | 19.2 | 54 | Average 3,000 | 1,000 | — | — |
| **CAMI II Toy Mouse Gut (short read)** | 26.7 | 64 | 2 × 150 | — | 270 | 20 |
| **Simulated Dataset** | 18.4 | 480 | 2 × 100 or 2 × 150 | — | 200, 300, 400 | 10, 25, 50 |

## Database Construction

The databases used in this experiment were previously built in the **[Build experiment](../build/README.md)**, utilizing the **completeONE** and **complete** databases for classification.

# 3. Environment Setup

## Software Requirements

### Chimera Version

The **Chimera** version used in this experiment is **v1.6.0**. It is recommended to use this version or a higher version to ensure compatibility and optimal performance.

### Other Dependencies

In addition to **Chimera**, the following software and libraries are required for this experiment:

- **multitax**: The version used in the experiment is **v1.3.1**, available at [multitax GitHub](https://github.com/pirovc/multitax).
- **numpy**: Used for numerical computations, the version used is **v1.24.4**.
- **Biopython**: Used for bioinformatics data processing, the latest version is recommended.

## Comparison Software

This experiment compares **Chimera** with the following classification tools:

- **Kraken2 v2.1.3**: [GitHub - Kraken2](https://github.com/DerrickWood/kraken2)
- **Bracken v2.9**: [GitHub - Bracken](https://github.com/jenniferlu717/Bracken)
- **Ganon v2.1.0**: [GitHub - Ganon](https://github.com/pirovc/ganon)
- **Taxor v0.1.3**: [GitHub - Taxor](https://github.com/JensUweUlrich/Taxor)

These tools are all used for metagenomic classification and will be compared with **Chimera** on various experimental metrics to assess their performance in terms of accuracy, efficiency, and resource consumption.

## Hardware Requirements

To ensure smooth execution of the experiment, the following hardware configurations are recommended:

- **CPU**: At least **8 cores** or more for optimal performance. More cores will help speed up processing, especially when handling large datasets and databases.
- **Memory**: At least **128 GB RAM**, especially due to the large size of the **Kraken2** database. Other classification tools and datasets also require significant memory.
- **Storage**: At least **500 GB of available storage** to accommodate the large databases and test datasets used by the various tools. Sufficient storage is crucial for handling the large size of the classification databases.

Depending on the dataset size and complexity of the experiment, more powerful hardware may be necessary.

# 4. Experimental Method

## Classification Experiment Process

In this experiment, the classification tasks are executed using **[runClassifier.py](../runClassifier.py)**. To handle multiple datasets more efficiently, we use **[generate_and_run.py](../generate_and_run.py)** to call **runClassifier.py**, simplifying the experiment workflow.

### Execution Command

Each software tool uses the same **70% classification threshold** and **32 threads**, with all other parameters set to their default values. The following are the execution commands for each classification tool:

```bash
# Chimera
conda activate chimera
/usr/bin/time -v chimera classify -i {seq_path} -d {db_path}DB.imcf -t {threads} -s {threshold} -o {output}

# Ganon
conda activate ganon
/usr/bin/time -v ganon classify -d {db_path} -s {seq_path} -t {threads} -c {threshold} -o {output} --verbose --output-all --output-one

# Ganon2
conda activate ganon
/usr/bin/time -v ganon classify -d {db_path} -s {seq_path} -t {threads} -c {threshold} -o {output} --verbose --output-all --output-one

# Kraken2
conda activate kraken2
/usr/bin/time -v kraken2 --db {db_path} --threads {threads} --confidence {threshold} --output {output} {extra} {seq_path} --report {report}

# Taxor
conda activate taxor
/usr/bin/time -v taxor search --index-file {db_path}.hixf --query-file {seq_path} --output-file {output} --percentage {threshold} --threads {threads}

# Bracken
conda activate kraken2
/usr/bin/time -v bracken -d {db_path} -i {report} -o {bracken_output}
```
### Important Parameters

- **Classification Threshold**: All tools use a 70% classification threshold.
- **Threads**: The default is set to 32 threads, as Taxor supports a maximum of 32 threads, and the same is applied to all tools in this experiment.