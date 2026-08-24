# Quality Control Images

This page contains the quality-control plots generated during the Seurat analysis of the Kramnik mouse saline Day 17 sample.

## Pre-filtering QC

### Violin Plots

![Pre-filtering QC violin plots](pre_filter_violin_plot.png)

These violin plots show the distributions of `nFeature_RNA`, `nCount_RNA`, and `percent.mt` before QC filtering.

### Feature Scatter Plots

![Pre-filtering QC scatter plots](pre_filter_feature_scatter.png)

These scatter plots show the relationships between:

- `nCount_RNA` and `percent.mt`
- `nCount_RNA` and `nFeature_RNA`

## Post-filtering QC

### Violin Plots

![Post-filtering QC violin plots](post_filter_violin_plot.png)

These violin plots show the QC metric distributions after applying the selected filtering thresholds.

### Feature Scatter Plots

![Post-filtering QC scatter plots](post_filter_feature_scatter.png)

These scatter plots show the relationships between the QC metrics after filtering.

## QC Summary

| QC metric | Threshold |
|---|---|
| **nFeature_RNA** | > 200 |
| **nFeature_RNA** | < 10,000 |
| **percent.mt** | < 5% |

For this tutorial run:

| QC result | Number of cells |
|---|---:|
| **Cells before QC** | 5,000 |
| **Cells after QC** | 5,000 |
| **Cells removed** | 0 |

No cells were removed using the selected QC thresholds in this 5,000-cell tutorial subset.
