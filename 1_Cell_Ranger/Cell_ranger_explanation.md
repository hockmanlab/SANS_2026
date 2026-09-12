# 1.1 Cell Ranger Count

After sequencing, you will receive **FASTQ files** containing the raw sequencing reads from your samples. These files need to be processed before we can perform downstream single-cell RNA-seq analysis in Seurat.

For data generated using 10x Genomics platforms, we can use **Cell Ranger**, a command-line software suite developed by 10x Genomics. In this tutorial, we will use the `cellranger count` pipeline.

Cell Ranger takes the FASTQ files and a reference transcriptome and processes the sequencing reads to generate **gene expression counts for individual cells**. Briefly, it identifies cell barcodes and UMIs, aligns the sequencing reads to the reference genome/transcriptome, and assigns the reads to genes.

The resulting gene-by-cell expression matrix is one of the main inputs we will use to create our Seurat object.

---

## Running Cell Ranger

For our data, Cell Ranger was run using the following Bash script:

```bash
#!/bin/bash

cellranger count --id=cell_ranger_d17_saline_KJ1 \
   --fastqs=/scratch/jmlkha003/kramnik/kramnik_samples_2025/d17_saline_KJ1/KJ1 \
   --sample=KJ1-SCI7T008-SCI5T008_2373HCLT4 \
   --transcriptome=/scratch/jmlkha003/kramnik_samples/alignment_files/refdata-gex-GRCm39-2024-A \
   --include-introns=true
```

The main arguments in this command are:

* `--id=run_count_MJ5`  
  Specifies the name of the Cell Ranger run and the directory in which the output will be stored.

* `--fastqs=/scratch/jmlkha003/MJ5`  
  Specifies the location of the FASTQ files that will be processed.

* `--transcriptome ../refdata-gex-mm10-2020-A`  
  Specifies the **10x-compatible reference transcriptome** used to align and annotate the sequencing reads. Here, the reference is the mouse **mm10** reference.

* `--sample=...`  
  Specifies the FASTQ sample prefixes that should be processed. Multiple sample prefixes can be provided as a comma-separated list.

* `--include-introns=true`  
  Tells Cell Ranger to include reads that map to **intronic regions** when calculating gene expression counts.

---

## Understanding File and Directory Paths

When running Cell Ranger, you will need to provide the locations of your FASTQ files and reference transcriptome. These locations are written as **directory paths**.

For example:

```bash
--fastqs=/scratch/jmlkha003/MJ5
```

This tells Cell Ranger where to find the FASTQ files.

The path:

```text
/scratch/jmlkha003/MJ5
```

is an example of a directory path on the **UCT HPC system**. The `/scratch/` directory is a storage location on the HPC system, `jmlkha003` is the user's HPC directory, and `MJ5` is the folder containing the FASTQ files.

For more information about the **UCT High Performance Computing (HPC) system**, see the [UCT HPC Cluster documentation](https://ucthpc.uct.ac.za/index.php/hpc-cluster/).

Your own directory path will most likely be different.

For example, if your FASTQ files are stored in another folder, you would replace:

```bash
--fastqs=/scratch/jmlkha003/MJ5
```

with the path to your own FASTQ directory:

```bash
--fastqs=/path/to/your/FASTQ_directory
```

Similarly, the reference transcriptome path may need to be changed depending on where the reference genome is stored on your computer or HPC system.

### Important

You should **not copy the example paths exactly** unless the files are actually located in those directories.

The general structure is:

```text
--fastqs=<directory containing your FASTQ files>
--transcriptome=<directory containing your reference transcriptome>
```

For example:

```bash
--fastqs=/path/to/my_sample \
--transcriptome=/path/to/my_reference
```

The `\` at the end of a line tells Bash that the command continues onto the next line. It is not part of the directory path.

---

## How to Find a Directory Path

If you are unsure where a file or folder is located, you can find its directory path using your operating system.

### macOS

On a Mac, you can find the path to a folder using **Finder**.

1. Open **Finder**.
2. Navigate to the folder containing your FASTQ files or reference.
3. Right-click the folder.
4. Select **Get Info** to view information about the folder.

You can also hold the **Option (⌥)** key while right-clicking the folder and select **Copy "... " as Pathname**. This copies the full path to your clipboard.

You can then paste the path into your R or Bash code.

For example:

```text
/Users/yourname/Documents/SANS_2026/data/MJ5
```

This could then be used as:

```bash
--fastqs=/Users/yourname/Documents/SANS_2026/data/MJ5
```

Another option is to open **Terminal**, type `cd ` (including the space), and then drag the folder from Finder into the Terminal window. macOS will automatically insert the folder's path.

### Windows

On Windows, you can find the path to a folder using **File Explorer**.

1. Open **File Explorer**.
2. Navigate to the folder containing your FASTQ files or reference.
3. Click the **address bar** at the top of the window.
4. The folder path will be displayed.
5. You can copy the path from the address bar.

For example:

```text
C:\Users\yourname\Documents\SANS_2026\data\MJ5
```

When using paths in R, it is often easier to use forward slashes:

```r
"C:/Users/yourname/Documents/SANS_2026/data/MJ5"
```

If you are running Bash through **WSL (Windows Subsystem for Linux)**, the path format is different. Windows drives are accessed through `/mnt/`.

For example:

```text
C:\Users\yourname\Documents\SANS_2026\data\MJ5
```

would typically be accessed in WSL as:

```text
/mnt/c/Users/yourname/Documents/SANS_2026/data/MJ5
```

> **Note:** The directory path depends on where you saved your files and which operating system or computing environment you are using.

---

## Why Include Intronic Reads?

In a single-cell RNA-seq experiment, not all captured RNA molecules correspond to fully spliced mature mRNA. Some reads originate from **unspliced transcripts** and therefore map to intronic regions.

Including these reads can increase the amount of usable data and improve the **sensitivity** of the analysis, resulting in higher UMI counts and more genes detected per cell.

This is particularly relevant for **single-nucleus RNA-seq**, where unspliced transcripts and intronic reads are often abundant.

For Cell Ranger **v7.0 and later**, intronic reads are included by default for whole-transcriptome Gene Expression data. We explicitly set `--include-introns=true` here to make the analysis setting clear and reproducible.

---

# 1.2 Cell Ranger Count Outputs

After Cell Ranger finishes, the `outs/` directory contains several files and folders. These provide processed expression matrices, quality-control metrics, and files that can be used for further exploration of the data.

A typical output directory will contain:

```text
run_count_MJ5/
└── outs/
    ├── analysis/
    ├── molecule_info.h5
    ├── cloupe.cloupe
    ├── raw_feature_bc_matrix/
    ├── filtered_feature_bc_matrix/
    ├── raw_feature_bc_matrix.h5
    ├── filtered_feature_bc_matrix.h5
    ├── web_summary.html
    └── metrics_summary.csv
```

The exact files and folders can vary depending on the Cell Ranger version and the options used when running the analysis.

---

## Cell Ranger Web Summary

One of the most useful Cell Ranger outputs is:

```text
web_summary.html
```

This is an interactive HTML report that summarises the Cell Ranger run and provides important **quality-control metrics**.

The HTML file contains information such as:

* sequencing and mapping statistics
* estimated number of cells
* sequencing saturation
* median genes detected per cell
* median UMIs detected per cell
* reads mapped to the genome
* reads mapped to transcriptome features
* other quality-control metrics

### Viewing the Example Cell Ranger Web Summary

The Cell Ranger output file for this tutorial is available from the Hockman Lab GitHub repository.

**The HTML file is too large for GitHub to display directly in the browser.**

You can download the `saline_day17_web_summary.html` file using the link below and then open it in your web browser to view the interactive Cell Ranger report:

[Download the `saline_day17_web_summary.html` file](https://github.com/hockmanlab/SANS_2026/blob/main/1_Cell_Ranger/saline_day17_web_summary.html)

After downloading the file:

1. Locate the file in your **Downloads** folder.
2. Double-click the HTML file.
3. The file should open in your default web browser.
4. You can then explore the Cell Ranger quality-control metrics and summary information.

> **Note:** GitHub may not display this large HTML file directly. Downloading the file and opening it locally in a web browser allows the interactive report to be viewed.

---

## Main Cell Ranger Output Files

The main output files are:

* **`analysis/`** — contains additional Cell Ranger analysis results and outputs, including clustering and dimensional-reduction results.

* **`molecule_info.h5`** — contains detailed information about detected molecules, including cell barcodes, UMIs, and associated features. This file can be used for more detailed downstream analysis.

* **`cloupe.cloupe`** — a file that can be opened in **Loupe Browser** to interactively explore the single-cell dataset, including cell clusters and gene expression.

* **`raw_feature_bc_matrix/`** — contains the unfiltered gene expression matrix. It includes barcodes and their associated expression counts before Cell Ranger applies cell-calling/filtering.

* **`filtered_feature_bc_matrix/`** — contains the filtered gene expression matrix for the cells identified by Cell Ranger.

* **`raw_feature_bc_matrix.h5`** — contains the raw expression matrix stored as a single HDF5 file.

* **`filtered_feature_bc_matrix.h5`** — contains the filtered expression matrix stored as a single HDF5 file. It contains the same expression information as the three files inside the `filtered_feature_bc_matrix/` folder, but stores the information in a single file.

* **`web_summary.html`** — an interactive HTML report summarising the Cell Ranger run and providing important quality-control metrics.

* **`metrics_summary.csv`** — contains key Cell Ranger summary metrics in a tabular format, which can be useful for checking the overall quality of the sequencing and cell-calling results.

---

## The Filtered Expression Matrix

For this tutorial, we will use the **`filtered_feature_bc_matrix.h5` file** as the main input for loading the Cell Ranger data into Seurat.

```text
filtered_feature_bc_matrix.h5
```

The `.h5` file contains the **filtered gene expression matrix** in a single HDF5 file. It contains the same expression information as the three files stored inside the `filtered_feature_bc_matrix/` folder, but stores the information together in one file.

Using the `.h5` file is convenient because you only need to specify the path to a single file when loading the data.

For example:

```r
library(Seurat)

DATA_FILE <- "/path/to/filtered_feature_bc_matrix.h5"

query_10x <- Read10X_h5(
  DATA_FILE
)
```

The path to the `.h5` file will need to be changed to the location of the file on your computer or HPC system.

### Alternatively: Read the Matrix Folder

If you want, you can also read the expression data directly from the `filtered_feature_bc_matrix/` folder.

The folder contains three files:

```text
filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

The three files work together to represent the **gene-by-cell expression matrix**.

* **`barcodes.tsv.gz`** — contains the cell barcodes, with each barcode representing an individual cell or nucleus.

* **`features.tsv.gz`** — contains information about the genes/features detected in the experiment, including gene IDs and gene names.

* **`matrix.mtx.gz`** — contains the actual expression counts, representing the number of detected UMIs for each gene in each cell or nucleus.

The folder can be read into Seurat using:

```r
library(Seurat)

DATA_DIR <- "/path/to/filtered_feature_bc_matrix"

query_10x <- Read10X(
  data.dir = DATA_DIR
)
```

Both approaches provide the filtered gene expression matrix needed for downstream Seurat analysis. For this tutorial, we will use the **`.h5` file as the main option**, but the `filtered_feature_bc_matrix/` folder can also be used if preferred.

---

> **Note:** The exact output files can vary depending on the Cell Ranger version and options used. Always check the documentation for the version of Cell Ranger being used.

For the official Cell Ranger command-line documentation, see the [10x Genomics Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-command-line-arguments).
