# Build Directory README

This README explains how the experiments in the `build` directory are conducted, detailing the process for running build-related benchmarks for each software included in ChimeraBenchmark. Additionally, it introduces the datasets used for database construction to ensure clarity and reproducibility.

---

## 1. Datasets Used for Database Construction

The following datasets were used for the experiments, all downloaded using **[Genome Updater](https://github.com/pirovc/genome_updater)**

| Dataset Name              | Size (GB) | Total Sequences | Base Pairs (BP)   | Assemblies |
|---------------------------|-----------|-----------------|-------------------|------------|
| Archaea (2024.10.10)      | 1.6       | 1138            | 1,725,129,521     | 601        |
| CompleteONE (2024.9.26)   | 50.6      | 36,811          | 52,322,440,925    | 25,091     |
| Complete (2024.10.7)      | 179.4     | 119,259         | 189,361,721,135   | 58,672     |
| Refseq (2024.10.16)       | 3248.3    | 51,704,208      | 3,430,607,682,023 | 403,253    |

These datasets were selected to evaluate each tool's performance on a range of sizes and complexities, from small datasets like `Archaea` to the massive `Refseq`.

---

## 2. [Chimera](https://github.com/LoadStar822/Chimera)
**Command**:
```bash
chimera build -i target.tsv -o ChimeraDB -m normal -t 32
```

**Description**:
- **Input**: `target.tsv` specifies the sequences and taxonomic identifiers.
- **Mode**: The `-m normal` mode constructs a 16-bit interleaved cuckoo filter for higher accuracy. The `fast` mode can also be used for quicker builds but with reduced accuracy.
- **Threads**: `-t 32` utilizes 32 threads for faster database construction.
- The resulting database is saved as `ChimeraDB`.

---

## 3. [Kraken2](https://github.com/DerrickWood/kraken2) and [Bracken](https://github.com/jenniferlu717/Bracken)
This section provides step-by-step instructions for constructing databases using **Kraken2** and **Bracken**, specifically tailored for the `Archaea` dataset. The same approach can be adapted for other datasets by changing the database name.

---

### 1. Kraken2 Database Construction

**Step-by-Step Commands**:

1. **Download the Taxonomy Dump**:
   ```bash
   kraken2-build --download-taxonomy --db Archaea
   ```

2. **Copy Sequence Files to the Kraken2 Directory**:
   ```bash
   find ~/tianqinzhong/project/chimera/archaea_files -type f -name '*.fna.gz' -exec cp {} ~/tianqinzhong/project/kraken2/Archaea/files/ \;
   ```

3. **Decompress the Files**:
   ```bash
   find . -maxdepth 1 -name "*.gz" -print0 | xargs -0 -P32 -n1 gunzip
   ```

4. **Add Files to the Kraken2 Library**:
   ```bash
   find ~/tianqinzhong/project/kraken2/Archaea/files -name "*.fna" -print0 | xargs -0 -n1 -I {} kraken2-build --add-to-library {} --db Archaea
   ```

5. **Build the Kraken2 Database**:
   ```bash
   /usr/bin/time -v kraken2-build --build --db Archaea --threads 32
   ```

6. **Repeat for Other Databases**:
   To build other databases like `CompleteONE`, `Complete`, or `Refseq`, replace `Archaea` with the respective database name in the above commands.

---

### 2. Bracken Database Construction

**Step-by-Step Commands**:

1. **Set Shell Stack Size** (to avoid parameter list too long errors):
   ```bash
   ulimit -s 1024000
   ```

2. **Build the Bracken Database**:
   ```bash
   /usr/bin/time -v bracken-build -d Archaea -t 32 -k 35 -l 100
   ```

**Explanation of Parameters**:
- `-d Archaea`: Specifies the Kraken2 database to use as input.
- `-t 32`: Sets the number of threads for parallel processing.
- `-k 35`: Sets the k-mer length for Bracken.
- `-l 100`: Specifies the read length for the Bracken database.

---

### Notes
- The commands provided are for the `Archaea` dataset. For other datasets, replace `Archaea` with the desired database name.
- Ensure that all necessary sequence files are available and correctly organized in the `files` directory.
- The use of `ulimit` ensures smooth execution when dealing with large numbers of input files.

By following these steps, you can construct Kraken2 and Bracken databases for various datasets and efficiently benchmark their performance.

---

## 4. [Ganon and Ganon2](https://github.com/pirovc/ganon)
This section provides instructions for constructing databases using **Ganon** and **Ganon2**. The primary difference between the two is that **Ganon** uses the `-v ibf` parameter to construct interleaved Bloom filters, whereas **Ganon2** omits this parameter. The following commands can be adapted for either tool by including or excluding `-v ibf`.

---

### 1. Database Construction Commands

#### **CompleteONE Database**
Constructed using RefSeq data with selected organism groups:
```bash
/usr/bin/time -v ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 32 --complete-genomes --genome-updater "-A 'species:1'" --db-prefix completeONE -v ibf
```

#### **Archaea Database**
Constructed using RefSeq data, limited to archaea:
```bash
/usr/bin/time -v ganon build --source refseq --organism-group archaea --threads 32 --complete-genomes --db-prefix Archaea -v ibf
```

#### **RefSeq Database**
Comprehensive database using RefSeq with all organism groups:
```bash
/usr/bin/time -v ganon build --source refseq --threads 32 --db-prefix refseq --organism-group archaea bacteria fungi human invertebrate metagenomes other plant protozoa vertebrate_mammalian vertebrate_other viral -v ibf
```

#### **Complete Database**
Constructed using RefSeq data for complete genomes of selected organism groups:
```bash
/usr/bin/time -v ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 32 --complete-genomes --db-prefix complete -v ibf
```

---

### 2. Adjusting for Ganon2
To use **Ganon2**, simply remove the `-v ibf` parameter from the above commands. For example:

**Ganon2 Command for CompleteONE Database**:
```bash
/usr/bin/time -v ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 32 --complete-genomes --genome-updater "-A 'species:1'" --db-prefix completeONE
```

---

### 3. Explanation of Parameters
- `--source refseq`: Specifies that the database source is NCBI RefSeq.
- `--organism-group`: Defines the organism groups to include in the database.
- `--threads 32`: Uses 32 threads for parallel processing.
- `--complete-genomes`: Limits the database to complete genomes (optional, based on the experiment).
- `--genome-updater "-A 'species:1'"`: Filters entries to include only species-level classifications.
- `--db-prefix`: Sets the prefix for the output database files.
- `-v ibf`: Constructs interleaved Bloom filters (used only in Ganon).

---

### 4. Notes
- The commands provided can be customized for different datasets by adjusting the `--organism-group` and `--db-prefix` parameters.
- Use **Ganon** for IBF-based filtering or **Ganon2** for standard Bloom filters by toggling the `-v ibf` option.
- Ensure sufficient system resources (e.g., memory and CPU) when running these commands, especially for large datasets like RefSeq.

By following these steps, you can construct Ganon and Ganon2 databases for various datasets and efficiently benchmark their performance.


---

## 5. [Taxor](https://github.com/JensUweUlrich/Taxor)
This section explains how to construct a database using **Taxor** with its default parameters. The instructions are tailored for the `CompleteONE` dataset, but the process can be adapted for other datasets by changing the names and paths accordingly.

---

### 1. Database Construction Commands

**Step-by-Step Instructions for CompleteONE**:

1. **Prepare the Taxonomy Dump**:
   Unpack the taxonomy dump into a dedicated folder:
   ```bash
   mkdir -p taxdump
   tar -zxvf taxdump.tar.gz -C taxdump
   ```

2. **Create a TSV File for Taxonomic Information**:
   Generate a `refseq_accessions_taxonomy.csv` file using `cut` and `taxonkit`:
   ```bash
   cut -f 1,7,20 assembly_summary.txt \
   | taxonkit lineage -i 2 -r -n -L --data-dir taxdump \
   | taxonkit reformat -I 2 -P -t --data-dir taxdump \
   | cut -f 1,2,3,4,6,7 > refseq_accessions_taxonomy.csv
   ```

3. **Build the Taxor Database**:
   Use the `taxor build` command to construct the database:
   ```bash
   /usr/bin/time -v taxor build --input-file completeONE.csv \
   --input-sequence-dir completeONE/2024-09-26_11-57-14/files \
   --output-filename completeONE.hixf \
   --threads 32 --kmer-size 22 --syncmer-size 12 --use-syncmer
   ```

**Explanation of Parameters**:
- `--input-file`: The path to the input CSV file (e.g., `completeONE.csv`).
- `--input-sequence-dir`: The directory containing all sequence files. Ensure all files are in a single directory without subfolders.
- `--output-filename`: The name of the output database file.
- `--threads 32`: The number of threads for parallel processing.
- `--kmer-size 22`: Specifies the k-mer size.
- `--syncmer-size 12`: Specifies the syncmer size.
- `--use-syncmer`: Enables the use of syncmers for indexing.

---

### 2. Notes and Adjustments for Other Datasets

- **Change File and Directory Names**:
  Replace `completeONE` with the appropriate dataset name (e.g., `Archaea`, `Complete`, `RefSeq`) in all commands.

- **Sequence Files in One Directory**:
  Ensure all downloaded sequence files are placed in a single directory without subfolders. This is required for `taxor build` to process the sequences correctly.

- **Using Genome Updater**:
  When downloading sequences with `genome-updater`, use the `-a` option to keep the current version of the taxonomy database in the output folder:
  ```bash
  genome-updater.sh -a ...
  ```

---

### 3. Example Adjustments for Other Datasets

- **Archaea Dataset**:
  ```bash
  /usr/bin/time -v taxor build --input-file archaea.csv \
  --input-sequence-dir archaea/files \
  --output-filename archaea.hixf \
  --threads 32 --kmer-size 22 --syncmer-size 12 --use-syncmer
  ```

- **Complete Dataset**:
  ```bash
  /usr/bin/time -v taxor build --input-file complete.csv \
  --input-sequence-dir complete/files \
  --output-filename complete.hixf \
  --threads 32 --kmer-size 22 --syncmer-size 12 --use-syncmer
  ```

- **RefSeq Dataset**:
  ```bash
  /usr/bin/time -v taxor build --input-file refseq.csv \
  --input-sequence-dir refseq/files \
  --output-filename refseq.hixf \
  --threads 32 --kmer-size 22 --syncmer-size 12 --use-syncmer
  ```

By following these steps and ensuring all sequence files are correctly organized, you can successfully construct Taxor databases for various datasets.

---

## Notes
- **Standardization**: All tools are benchmarked using the same sequence files and taxonomic assignments to ensure fair comparisons.
- **Performance Metrics**: The time taken and resources consumed during the database build process are recorded for each tool.
- **Reproducibility**: Ensure that the input files (`target.tsv`, `target.txt`, etc.) match the ones described above for consistent results.

This README serves as a guide for setting up and running build experiments with the various tools included in ChimeraBenchmark.
