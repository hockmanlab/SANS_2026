# Clustering and Annotation

In this section, we will identify groups of cells with similar gene expression profiles
and determine their biological identities using marker genes and reference-based annotation.

## Workshop files

The following files are available in this repository:

-   `clustering_annotation_BA.R` -- R script for hands-on session
-   `clustering_annotation_BA.Rmd` -- Markdown script for hands-on
-   `clustering_annotation.html` -- Completed tutorial with expected
    outputs

The data files required for the analysis are too large to host directly in this repository and **will be provided during the workshop**:

* `filtered_feature_bc_matrix.h5` – 10x Genomics gene-expression count matrix
* `seurat_preprocessed_SCT.rds` – preprocessed Seurat object for the clustering and annotation practical

## Instructions

1. Download the workshop folder from this repository.
2. The required `.h5` and `.rds` data files will be provided during the workshop.
3. Place the provided data files in the appropriate workshop folder as instructed.
4. Open `clustering_annotation.R` in RStudio.
5. Install the required packages using the installation section of the script.
6. Run the script sequentially from top to bottom.

If you are unable to run the analysis during the practical, open `clustering_annotation.html` to follow the workflow and view the expected results.
