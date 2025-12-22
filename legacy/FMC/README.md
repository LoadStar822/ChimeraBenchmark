# FairMin-Cap Directory README

## 1. Testing Classification with Different Maximum Hash Counts

### Experiment Objective

This experiment aims to test the classification performance of **Chimera** with different **maximum hash counts (`--max-hashes`)** settings, evaluating the impact of the **`--max-hashes`** parameter on classification tasks. Specifically, we will explore the effect of varying hash counts (from 1 million to 100 million, finally set to 0) on classification accuracy, efficiency, and memory consumption.

### Experiment Method

We will adjust the **`--max-hashes`** parameter in the **Chimera** tool from **1 million** to **100 million**, and finally set it to **0**, to build the **complete** database. The **CAMI II Marine long-read dataset** will be used as the test dataset for classification.

1. **Build Database**: Use different **`--max-hashes`** parameters to construct the **complete** database.
2. **Perform Classification**: Use the constructed database to classify the **CAMI II Marine long-read dataset**, with classification tasks executed through the **[runClassifier.py](../runClassifier.py)** script.

### Command Example

#### Database Construction

For example, to build the **complete** database with a specific **`--max-hashes`** value:

`/usr/bin/time -v chimera build -i complete/target.tsv -o completeDB -t 32 --max-hashes 1000000`

#### Classification

Classification tasks will be executed by calling the **[runClassifier.py](../runClassifier.py)** script, which will use the constructed database to perform the classification.

### Expected Output

The following results are expected from executing these commands:

- **Classification Accuracy**: Evaluation of classification accuracy, precision, recall, and other metrics for different **`--max-hashes`** settings.
- **Number of Truncated Species**: Comparison of the number of truncated species for different **`--max-hashes`** settings, analyzing the effect of database size on species capture.
- **Database Size**: Evaluation of the size of the **complete** database. As the **`--max-hashes`** parameter increases, the database size will increase accordingly.
- **Performance Metrics**: Comparison of database construction time, classification time, and memory consumption for different hash count settings.
- **Resource Consumption**: Observing the impact of different hash counts on system resources (e.g., memory, CPU usage).

## 2. Comparison of Classification Performance with 2 Million and 4 Million Hashes

### Experiment Objective

This experiment aims to compare the classification performance of **Chimera** using **2 million** and **4 million** maximum hash counts, specifically evaluating the performance of the **FairMin-Cap (FMC)** strategy on a larger and more complex **NCBI Refseq** database. By testing different **`--max-hashes`** settings on the **NCBI Refseq** database, we analyze the impact of the **FMC** strategy on classification tasks when handling a larger-scale dataset.

### Experiment Method

We will use **2 million** and **4 million** different **`--max-hashes`** settings to build the **NCBI Refseq** database, and use the **CAMI II Marine long-read dataset** as the test dataset for classification tasks.

1. **Build Database**: Construct the **NCBI Refseq** database with **2 million** and **4 million** hash settings.
2. **Perform Classification**: Use the built **NCBI Refseq** database to classify the **CAMI II Marine long-read dataset** and evaluate the performance of the **FMC strategy** on larger databases.

### Command Example

#### Database Construction

Use the following commands to build the **NCBI Refseq** database, setting **`--max-hashes`** to **2 million** and **4 million**:

`/usr/bin/time -v chimera build -i NCBIRefseq/target.tsv -o refseqDB_2M -t 32 --max-hashes 2000000`

`/usr/bin/time -v chimera build -i NCBIRefseq/target.tsv -o refseqDB_4M -t 32 --max-hashes 4000000`

#### Classification

Classification tasks will be executed by calling the **[runClassifier.py](../runClassifier.py)** script.

### Expected Output

By executing these commands, we expect the following results:

- **Classification Accuracy**: Evaluate the classification accuracy, precision, recall, and other metrics for different **`--max-hashes`** settings.
- **Performance Evaluation**: Compare database construction time, classification time, and memory consumption for different hash count settings, and assess the performance of the **FMC strategy** on large-scale databases.
- **Resource Consumption**: Observe the impact of different hash counts on system resources (e.g., memory, CPU), and evaluate how the **FMC strategy** affects resource consumption.


## 3. Modifying Ganon Source Code to Apply FMC Strategy

### Experiment Objective

This experiment modifies the **Ganon** source code in **GanonBuild.cpp** to apply the **FairMin-Cap (FMC)** strategy and compares the performance of the modified **Ganon** and the original **Ganon** in classification tasks. The goal is to assess the effectiveness of the **FMC strategy** in **Ganon**, especially for large-scale databases.

### Experiment Method

We modified the **GanonBuild.cpp** source code to integrate the **FMC strategy**. The modified source file is available in the current directory, and can be accessed via the following link: **[GanonBuild.cpp](./GanonBuild.cpp)**. We will use the **CAMI II Marine long-read dataset** as the test dataset, perform classification using both the modified and the original **Ganon**, and compare their performance.

1. **Modify the Source Code**: Modify the **GanonBuild.cpp** source code in **Ganon** to apply the **FMC strategy**. The modified file is available in the current directory.
2. **Build Database**: Use the modified **Ganon** to build the database. 
3. **Classification**: The classification task will be performed using the **[runClassifier.py](../runClassifier.py)** script, as done in previous experiments.

### Command Example

#### Database Construction

Use the following command to build the **Ganon** database:

`/usr/bin/time -v ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 32 --complete-genomes --db-prefix complete -v ibf`

#### Classification

The classification task will be executed by calling the **[runClassifier.py](../runClassifier.py)** script.

### Expected Output

By executing these commands, we expect the following results:

- **Classification Accuracy**: Evaluate the difference in classification accuracy, precision, recall, and other metrics between the modified **Ganon** and the original **Ganon**.
- **Performance Evaluation**: Compare the classification time, memory consumption, and other performance metrics between the modified **Ganon** and the original **Ganon**.
- **Resource Consumption**: Observe the impact of the **FMC strategy** on the efficiency and resource consumption of **Ganon**.

## 4. Random Selection of 54 Truncated Species for Comparison

### Experiment Objective

This experiment aims to randomly select **54 truncated species** that were significantly impacted under the **2 million hash setting** and evaluate the classification performance using the **FMC strategy**. Specifically, we want to test whether the **FMC strategy** has any negative effects on these species.

### Experiment Method

Under the **2 million hash setting**, a total of **1013 species** were affected. From these impacted species, we randomly selected **54 species** and created a test dataset (see file **[combined.fasta](./combined.fasta)**). This dataset is used to test whether the application of the **FMC strategy** negatively affects the classification of these truncated species.

1. **Select Species**: Randomly choose **54 species** from the **1013 affected species** and create the test dataset **combined.fasta**.
2. **Build Database**: Build the **complete database** with **2 million hash** and **`--max-hashes 0`** settings.
3. **Perform Classification**: Use the constructed databases to classify the **54 truncated species** and evaluate whether the **FMC strategy** negatively impacts these species.

### Command Example

#### Database Construction

Use the following command to build the **complete database**, with **`--max-hashes`** set to **2 million** and **0**:

`/usr/bin/time -v chimera build -i complete/target.tsv -o completeDB_2M -t 32 --max-hashes 2000000`

`/usr/bin/time -v chimera build -i complete/target.tsv -o completeDB_0 -t 32 --max-hashes 0`

#### Classification

The classification task will be executed using the **[runClassifier.py](../runClassifier.py)** script, and the **combined.fasta** dataset will be used for testing.

### Expected Output

By executing these commands, we expect the following results:

- **Classification Accuracy**: Evaluate the differences in classification accuracy, precision, recall, and other metrics between the **2 million hash complete database** and the **`--max-hashes 0` complete database**.
- **Performance Evaluation**: Compare database construction time, classification time, and memory consumption for different hash count settings.
- **Impact on Truncated Species**: Analyze how the **FMC strategy** affects the classification of the **54 truncated species** and assess whether it has any negative impact.
