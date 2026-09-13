# Cell Ranger Count

For the official Cell Ranger documentation, see the [10x Genomics Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct).

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

* `--id=cell_ranger_d17_saline_KJ1`  
  Specifies the name of the Cell Ranger run and the directory in which the output will be stored.

* `--fastqs=/scratch/jmlkha003/kramnik/kramnik_samples_2025/d17_saline_KJ1/KJ1`  
  Specifies the location of the FASTQ files that will be processed.

* `--transcriptome /scratch/jmlkha003/kramnik_samples/alignment_files/refdata-gex-GRCm39-2024-A`  
  Specifies the **10x-compatible reference transcriptome** used to align and annotate the sequencing reads. Here, the reference is the mouse **GRCm39** reference.

* `--sample=...`  
  Specifies the FASTQ sample prefixes that should be processed. Multiple sample prefixes can be provided as a comma-separated list.

* `--include-introns=true`  
  Tells Cell Ranger to include reads that map to **intronic regions** when calculating gene expression counts.

---

## Understanding File and Directory Paths

When running Cell Ranger, you will need to provide the locations of your FASTQ files and reference transcriptome. These locations are written as **directory paths**.

For example:

```bash
--fastqs=/scratch/jmlkha003/kramnik/kramnik_samples_2025/d17_saline_KJ1/KJ1
```

This tells Cell Ranger where to find the FASTQ files.

The path:

```text
/scratch/jmlkha003/kramnik
```

is an example of a directory path on the **UCT HPC system**. The `/scratch/` directory is a storage location on the HPC system, `jmlkha003` is the user's HPC directory, and `kramnik` is the folder containing the FASTQ files.

For more information about the **UCT High Performance Computing (HPC) system**, see the [UCT HPC Cluster documentation](https://ucthpc.uct.ac.za/index.php/hpc-cluster/).

Your own directory path will most likely be different.

For example, if your FASTQ files are stored in another folder, you would replace:

```bash
--fastqs=/scratch/jmlkha003/kramnik
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

You can then paste the path into your R or Bash code.

For example:

```text
/Users/yourname/Documents/SANS_2026/data/MJ5
```

This could then be used as:

```bash
--fastqs=/Users/yourname/Documents/SANS_2026/data/MJ5
```

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

> **Note:** The directory path depends on where you saved your files and which operating system or computing environment you are using.

---

## Why Include Intronic Reads?

In a single-cell RNA-seq experiment, not all captured RNA molecules correspond to fully spliced mature mRNA. Some reads originate from **unspliced transcripts** and therefore map to intronic regions.

Including these reads can increase the amount of usable data and improve the **sensitivity** of the analysis, resulting in higher UMI counts and more genes detected per cell.

This is particularly relevant for **single-nucleus RNA-seq**, where unspliced transcripts and intronic reads are often abundant.

For Cell Ranger **v7.0 and later**, intronic reads are included by default for whole-transcriptome Gene Expression data. We explicitly set `--include-introns=true` here to make the analysis setting clear and reproducible.

---

# Cell Ranger Count Outputs

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

## Cell Ranger Web Summary

One of the most useful Cell Ranger outputs is:

```text
web_summary.html
```

This is an interactive HTML report that summarises the Cell Ranger run and provides important **quality-control metrics**.

Some of the main metrics shown include:

- **Estimated number of cells:** The number of cell-associated barcodes identified by Cell Ranger.
- **Mean reads per cell:** The average number of sequencing reads associated with each identified cell.
- **Median genes per cell:** The median number of genes detected in each cell.
- **Median UMI counts per cell:** The median number of unique RNA molecules detected in each cell.
- **Sequencing saturation:** Indicates how much of the available library has been sampled by sequencing.
- **Fraction reads in cells:** The proportion of sequencing reads associated with identified cells rather than background droplets.
- **Reads mapped to genome:** The proportion of reads that align to the reference genome.
- **Reads mapped confidently to transcriptome:** The proportion of reads that can be confidently assigned to annotated transcripts.

These metrics provide an initial assessment of whether the library was successfully sequenced and whether the resulting data are suitable for downstream single-cell analysis.

### Cell Barcodes and UMIs

Two important concepts when interpreting the Cell Ranger metrics are **cell barcodes** and **UMIs**.

**Cell barcodes** are short DNA sequences used to identify which cell a sequencing read came from. Each droplet receives a barcode, allowing sequencing reads to be assigned back to individual cells.

**UMIs (Unique Molecular Identifiers)** are short sequences attached to individual RNA molecules. They allow Cell Ranger to distinguish original RNA molecules from duplicate copies produced during PCR amplification.

In simple terms:

- **Cell barcode → Which cell did the RNA come from?**
- **UMI → Which original RNA molecule is it?**
- **Gene → What RNA was detected?**

Therefore, metrics such as **UMIs per cell** provide an estimate of the number of unique RNA molecules detected from each cell.

## Knee Plot

The **knee plot** shows the distribution of UMI counts across the ranked cell barcodes.

- The **x-axis** shows cell barcodes ranked from the barcode with the highest number of UMIs to the barcode with the lowest number of UMIs.
- The **y-axis** shows the number of UMIs associated with each barcode.
- Barcodes with high UMI counts are more likely to represent **real cells**.
- Barcodes with very low UMI counts are more likely to represent **empty droplets or background RNA**.
- The **knee** is the region where the curve drops sharply and helps distinguish cell-associated barcodes from background.

> **In simple terms:** The knee plot helps show where the real cells end and the background/empty droplets begin.

Cell Ranger uses the barcode and UMI information, together with its cell-calling algorithm, to determine which barcodes are likely to represent cells.

### Viewing the Example Cell Ranger Web Summary

The Cell Ranger output file for this tutorial is available from the Hockman Lab GitHub repository.

**The HTML file is too large for GitHub to display directly in the browser.**

[Open the `saline_day17_web_summary.html` file on GitHub](https://github.com/hockmanlab/SANS_2026/blob/main/1_Cell_Ranger/saline_day17_web_summary.html)

Once the file opens on GitHub, click the **downward arrow (⬇) download button on the right-hand side** to download the HTML file.

After downloading:

1. Locate the `saline_day17_web_summary.html` file in your Downloads folder.
2. Double-click the file to open it in your web browser.
3. You can then explore the interactive Cell Ranger report and its quality-control metrics.

> **Note:** GitHub may not display the full interactive HTML report directly because of the file size. Downloading the file and opening it locally in a web browser allows you to view the report properly.

In the next section, we will read in the  **`filtered_feature_bc_matrix.h5`** and **`filtered_feature_bc_matrix/`** 

> **Note:** The exact output files can vary depending on the Cell Ranger version and options used. Always check the documentation for the version of Cell Ranger being used.

For the official Cell Ranger command-line documentation, see the [10x Genomics Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-command-line-arguments).
