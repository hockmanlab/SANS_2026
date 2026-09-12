Clustering and Annotation
================
Bismark Appiah, PhD.

In this section, we will identify groups of cells with similar gene
expression profiles and determine their biological identities using
marker genes and reference-based annotation.

The focus is on understanding the key concepts behind clustering and
annotation rather than simply running the commands.

### Pre-processed Seurat object

The pre-processed Seurat object can be downloaded here:

[Download the pre-processed Seurat
object](https://github.com/hockmanlab/SANS_2026/releases/download/SCTransform_preprocessed_data/seurat_preprocessed_SCT.rds)

## Workflow

The analysis follows this workflow:

# Clustering and Annotation

The analysis follows this workflow:

```text
Pre-processed Seurat object
        ↓
PCA
        ↓
Select principal PCs
        ↓
Find neighbours
        ↓
Clustering
        ↓
Compare resolutions
        ↓
UMAP
        ↓
Identify markers
        ↓
Manual annotation
        ↓
Reference annotation (Azimuth)
