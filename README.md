# ChimeraBenchmark

### 1. Project Overview

**ChimeraBenchmark** is a dedicated benchmarking project developed to support and evaluate the performance of **Chimera**, a high-speed metagenomic classification tool created by Qinzhong Tian. This project is designed to streamline the process of running experiments, analyzing results, and visualizing data, making it easier for researchers to understand the capabilities and efficiency of Chimera in comparison to other metagenomic tools.

**Key Features**:
- **Automated Software Pipeline**: ChimeraBenchmark can automatically run a series of popular metagenomic classification tools, including **Kraken2**, **Ganon**, **Ganon2**, **Taxor**, **Bracken**, and **Chimera**, given the necessary databases and test files. This automation simplifies the benchmarking process, enabling efficient and consistent testing across various tools.
- **Accuracy and Abundance Metrics**: If a ground truth file is provided, ChimeraBenchmark can automatically calculate the classification accuracy and compute the **L1 Distance** for abundance estimations, allowing for in-depth performance analysis.
- **Data Analysis and Visualization**: The project includes code for processing experiment outputs and generating visualizations, such as charts and graphs, to clearly present classification results and performance statistics.
- **Reproducibility**: ChimeraBenchmark ensures that all experiments can be reproduced easily, enabling consistent performance comparisons across different datasets and configurations.

**Use Cases**:
- **Performance Comparison**: Benchmark Chimera against other metagenomic classification tools to showcase its speed, accuracy, and memory efficiency.
- **Pipeline Customization**: Adapt the provided scripts for custom datasets or specific experimental setups to suit research needs.
- **Comprehensive Evaluation**: Use the automated pipeline to evaluate classification performance, measure accuracy, and assess the L1 Distance of abundance estimates with minimal manual intervention.
- **Publication Support**: Generate visual outputs suitable for presentations and publications, demonstrating Chimera’s performance metrics effectively.

With ChimeraBenchmark, researchers can confidently perform thorough analyses and share clear, visual results that highlight Chimera’s capabilities in large-scale metagenomic studies.

### 2. Quick Start

ChimeraBenchmark is designed to be simple and intuitive to use. By following the steps below, you can quickly run benchmarking tests for Chimera and other metagenomic classification tools.

**Basic Command Example**:
To run a basic benchmarking test, use the following command:

```bash
python -m runClassifier -s chimera kraken2 ganon -d /path/to/db_subpath -f /path/to/seq_file1.fasta /path/to/seq_file2.fasta -t 8 -o ./benchmark_results
```

**Explanation of Parameters**:
- `-s` or `--software`: Specifies the software to benchmark. You can include one or more tools from the following options: `chimera`, `ganon2`, `ganon`, `kraken2`, `taxor`, `bracken`.
- `-d` or `--db-subpath`: The subpath of the database to append to each software's base database path. Ensure the path is correct for the databases you are benchmarking.
- `-f` or `--seq-files`: The sequence file(s) to be processed. Multiple files can be specified, and both single and paired-end reads are supported.
- `-t` or `--threads`: The number of threads to use for parallel processing. The default is `4`, but higher numbers can be used for faster execution on multi-core systems.
- `-o` or `--output-dir`: The directory to store output results. If not specified, the default is `./benchmark_results`.
- `-sf` or `--standard-files`: (Optional) The path(s) to standard result files. If provided, ChimeraBenchmark will calculate classification accuracy and the L1 Distance for abundance analysis.

**Example Command with Standard Files**:
```bash
python -m runClassifier -s chimera ganon2 -d /path/to/db_subpath -f /path/to/seq_file.fasta -t 8 -o ./benchmark_results -sf /path/to/standard_result.txt
```

This command runs Chimera and Ganon2 using the specified database subpath and sequence file, outputs results to `./benchmark_results`, and calculates classification accuracy and L1 Distance using the provided standard result file.

### 3. Usage Guide

ChimeraBenchmark allows for flexible and comprehensive benchmarking with various configuration options. Below is a detailed breakdown of how to use the main parameters effectively:

**Running a Benchmark**:
1. **Select Software**: Use `-s` to choose the classification tools you want to benchmark. Multiple tools can be tested simultaneously.
2. **Specify Databases**: Ensure that each tool has access to its respective database by setting the `-d` parameter.
3. **Add Sequence Files**: Use `-f` to input one or more sequence files in FASTA format or other supported formats.
4. **Set Threads**: Increase `-t` to match the number of CPU cores available for faster execution.
5. **Output Directory**: Use `-o` to specify where results will be stored.
6. **Use Standard Files for Accuracy**: If you have ground truth data, include `-sf` to enable ChimeraBenchmark to calculate accuracy and the L1 Distance for abundance comparisons.

**Tips**:
- **Thread Management**: Adjust `-t` based on your machine’s capabilities for optimal performance.
- **File Management**: Ensure that the paths to sequence and standard result files are correct to prevent errors during execution.
- **Output Review**: Check the `./benchmark_results` directory (or the specified output path) for detailed logs and result files after running the benchmarks.

By following these steps and using the provided options, you can efficiently benchmark multiple metagenomic tools and analyze their performance.


### 4. Results Visualization

#### 4.1 Performance Metrics Visualization with `drawBoxLineDiagram.py`

ChimeraBenchmark includes functionality for visualizing performance metrics across different classification tools using the `drawBoxLineDiagram.py` script. This script, located in the `build` directory, automates the process of loading data from CSV files, processing it, and generating comprehensive performance plots.

**How It Works**:
- The script reads CSV files from a specified folder. Each CSV file should represent results for a specific tool and must include the following columns: `Dataset Name`, `Database`, `Taxonomic Rank`, `Total Samples`, `True Positives (TP)`, `False Positives (FP)`, `False Negatives (FN)`, `Accuracy`, `Precision`, `Recall`, and `F1 Score`.
- The script extracts the tool's name from the file name (e.g., `kraken2_results.csv` assumes "Kraken2").
- The script generates bar plots showing the average performance metrics (accuracy, precision, recall, F1 score) for each software across taxonomic ranks and databases.

**Usage Instructions**:
1. Run the script from the command line, specifying the folder containing the CSV result files:
   ```bash
   python drawBoxLineDiagram.py <folder_path>
   ```
   Replace `<folder_path>` with the path to the directory containing your CSV files.

2. Ensure that your CSV files are formatted correctly and contain the necessary columns.

**Example CSV Format**:
Each CSV file should have columns as shown below:

| Dataset Name         | Database     | Taxonomic Rank | Total Samples | True Positives (TP) | False Positives (FP) | False Negatives (FN) | Accuracy | Precision | Recall | F1 Score |
|----------------------|--------------|----------------|---------------|---------------------|----------------------|----------------------|----------|-----------|--------|----------|
| sample4_anonymous_gsa | completeONE | Superkingdom   | 148509        | 91308               | 10                   | 57191                | 0.6148   | 0.9999    | 0.6149 | 0.7615   |

**Visualization Details**:
- The script creates subplots for each combination of database and taxonomic rank.
- It displays bar plots with metrics represented by different colors.
- Each bar shows the average value for the given metric and software, with numerical labels for easy interpretation.
- The final output is saved as a high-resolution PDF file (`metrics_plot.pdf`) in the specified folder.

**Key Features**:
- **Comprehensive Overview**: Compares multiple tools across different taxonomic ranks and databases.
- **Customizable Appearance**: Utilizes Seaborn and Matplotlib for high-quality plots.
- **Ease of Use**: Simply place your CSV files in the folder, run the script, and get a visually rich PDF output.

**Output Example**:
- The generated PDF will contain a grid of plots organized by taxonomic rank and database, with each plot showing performance metrics for different tools.

**Requirements**:
- Ensure `pandas`, `seaborn`, and `matplotlib` are installed in your Python environment.
### 5. Support and Feedback
For any questions or support, feel free to reach out to us:
- **Website**: [MalabZ](http://lab.malab.cn/~cjt/MSA/)
- **Personal Homepage**: [Qinzhong Tian](https://loadstar822.github.io/)
- **Email**: tianqinzhong@qq.com
### 6. License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.