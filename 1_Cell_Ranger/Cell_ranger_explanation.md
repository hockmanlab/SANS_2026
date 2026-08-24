### Cell Ranger Count

After sequencing, you will receive **FASTQ files** containing the raw sequencing reads from your samples. These files need to be processed before we can perform downstream single-cell RNA-seq analysis in Seurat.

For data generated using 10x Genomics platforms, we can use **Cell Ranger**, a command-line software suite developed by 10x Genomics. In this tutorial, we will use the `cellranger count` pipeline. Cell Ranger takes the FASTQ files and a reference transcriptome and processes the sequencing reads to generate **gene expression counts for individual cells**. Briefly, it identifies cell barcodes and UMIs, aligns the sequencing reads to the reference genome/transcriptome, and assigns the reads to genes. The resulting gene-by-cell expression matrix is one of the main inputs we will use to create our Seurat object.

### Running Cell Ranger

For our data, Cell Ranger was run using the following Bash script:

```bash
#!/bin/bash

cellranger count --id=run_count_MJ5 \
--fastqs=/scratch/jmlkha003/MJ5 \
--transcriptome ../refdata-gex-mm10-2020-A \
--sample=MJ4-AK2831_H37WYDSX3,MJ4-AK2832_H37WYDSX3,MJ4-AK2833_H37WYDSX3,MJ4-AK2834_H37WYDSX3 \
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

### Why include intronic reads?

In a single-cell RNA-seq experiment, not all captured RNA molecules correspond to fully spliced mature mRNA. Some reads originate from **unspliced transcripts** and therefore map to intronic regions. 10x Genomics reports that intronic reads can represent a substantial proportion of reads in whole-transcriptome Gene Expression data.

Including these reads can increase the amount of usable data and improve the **sensitivity** of the analysis, resulting in higher UMI counts and more genes detected per cell. This is particularly relevant for **single-nucleus RNA-seq**, where unspliced transcripts and intronic reads are often abundant.

For Cell Ranger **v7.0 and later**, intronic reads are included by default for whole-transcriptome Gene Expression data. We explicitly set `--include-introns=true` here to make the analysis setting clear and reproducible.

### Cell Ranger Count Outputs

After Cell Ranger finishes, the `outs/` directory contains several files and folders. These provide the processed expression matrices, quality-control metrics, and files that can be used for further exploration of the data.

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

The main output files are:

* **`analysis/`** — contains additional Cell Ranger analysis results and outputs, including clustering and dimensionality-reduction results.
* **`molecule_info.h5`** — contains detailed information about the detected molecules, including their cell barcodes, UMIs, and associated features. This file can be used for more detailed downstream analysis.
* **`cloupe.cloupe`** — a file that can be opened in **Loupe Browser** to interactively explore the single-cell dataset, including cell clusters and gene expression.
* **`raw_feature_bc_matrix/`** — contains the unfiltered gene expression matrix. It includes barcodes and their associated expression counts before Cell Ranger applies its cell-calling/filtering.
* **`filtered_feature_bc_matrix/`** — contains the filtered gene expression matrix for the cells identified by Cell Ranger. This is one of the main outputs we will use for our Seurat analysis.
* **`raw_feature_bc_matrix.h5`** — the raw expression matrix stored as a single HDF5 file.
* **`filtered_feature_bc_matrix.h5`** — the filtered expression matrix stored as a single HDF5 file. It contains the same expression information as the three files inside `filtered_feature_bc_matrix/`, but in a single file.
* **`web_summary.html`** — an interactive HTML report summarising the Cell Ranger run and providing important quality-control metrics.
* **`metrics_summary.csv`** — contains key Cell Ranger summary metrics in a tabular format, which can be useful for checking the overall quality of the sequencing and cell-calling results.

### The filtered expression matrix

For this tutorial, we will mainly work with the **filtered** expression matrix:

```text
filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

The three files work together to represent the gene-by-cell expression matrix:

* **`barcodes.tsv.gz`** — contains the cell barcodes, with each barcode representing an individual cell.
* **`features.tsv.gz`** — contains information about the genes/features detected in the experiment, including gene IDs and gene names.
* **`matrix.mtx.gz`** — contains the actual expression counts, representing the number of detected UMIs for each gene in each cell.

Cell Ranger also provides:

```text
filtered_feature_bc_matrix.h5
```

This `.h5` file contains the same **filtered gene expression matrix** as the three files above, but stores the information together in a single HDF5 file. It is therefore a convenient alternative when loading the data into Seurat.

In the next section, we will show how to load the expression data into Seurat using **both formats**: the `filtered_feature_bc_matrix` folder and the `filtered_feature_bc_matrix.h5` file.


> **Note:** The exact output files can vary depending on the Cell Ranger version and options used. For example, newer Cell Ranger versions have additional requirements for BAM generation. Always check the documentation for the version being used.

For the official Cell Ranger command-line documentation, see [10x Genomics Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-command-line-arguments?utm_source=chatgpt.com).
