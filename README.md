# SANS_2026
# 🧬 Single-Nucleus RNA Sequencing Workshop

Welcome to an introduction to **single-nucleus RNA sequencing (snRNA-seq)** and the analysis of cellular heterogeneity at single-nucleus resolution.

This workshop combines an introduction to single-nucleus biology with a **hands-on exploration of a Kramnik mouse brain saline sample generated using 10x Genomics single-nucleus RNA sequencing technology**.

## 🔬 What is Single-Nucleus RNA Sequencing?

Single-nucleus RNA sequencing allows us to examine gene expression at the level of individual **cell nuclei**. It is particularly useful for studying complex tissues such as the **brain**, where isolating intact cells can be challenging.

Single-nucleus analysis allows us to explore:

- Cell-type composition
- Cellular heterogeneity
- Cell-specific gene expression
- Marker genes and cell-type identities and more.

## 🧪 Practical Dataset

During the hands-on session, we will explore a **Kramnik mouse brain saline sample** prepared using the **10x Genomics single-nucleus RNA-seq workflow**.

We will work through the main steps of a single-nucleus analysis workflow:

```text
10x Genomics snRNA-seq data
          ↓
Create Seurat object
          ↓
Quality control
          ↓
Normalisation
          ↓
Dimensionality reduction
          ↓
Clustering
          ↓
Identify marker genes
          ↓
Cell-type annotation
          ↓
Explore and interpret the data
```

## 💻 Hands-On Analysis

Participants will use **R and Seurat** to:

- Load 10x Genomics single-nucleus data
- Explore the dataset
- Perform quality control
- Normalise gene expression data
- Identify highly variable genes
- Explore nuclear populations using PCA and UMAP
- Identify clusters
- Find marker genes
- Investigate cell-type identities
- Visualise and interpret the data

## 🎯 Workshop Goal

By the end of the workshop, participants will have a basic practical introduction to how **10x Genomics single-nucleus RNA-seq data are processed and analysed**, from initial data exploration through to cell-type identification and biological interpretation.

The practical session will take participants **from 10x single-nucleus data to an exploration of the cellular landscape of the mouse brain**.

## 🧠 Dataset

| Information | Details |
|---|---|
| **Species** | Mouse |
| **Tissue** | Brain - cortex|
| **Condition** | Saline control |
| **Technology** | 10x Genomics Single-Nucleus RNA Sequencing (snRNA-seq) |
| **Analysis** | R / Seurat |

## 📋 Prerequisites

Before starting the workshop, please make sure that **R, RStudio, and all required R packages** are installed.

👉 **[Installation Guide](./2_Seurat/2.1_installation.md)**

Please complete the installation before starting the practical session.

## System requirements

The workshop was tested on macOS and Windows systems with 8 - 16 GB RAM.
A minimum of 8 GB RAM is recommended; 16 GB or more is preferable for smoother performance. Systems with less RAM may run more slowly or encounter memory-related errors.
