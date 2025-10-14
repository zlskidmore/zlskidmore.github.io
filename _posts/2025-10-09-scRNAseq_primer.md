---
layout: post
title:  "scRNAseq - Primer"
date:   2025-10-09 10:13:50 -0500
description: "A Primer on Single-Cell RNA-seq: Concepts, Steps, and Pitfalls"
tags: [scRNAseq, bioinformatics]
---

### Introduction
Welcome to a primer on single-cell RNA sequencing (scRNA-seq). Studies suggest that new information is retained more effectively when paired with a high-level overview, this overview act as a sort of cognitive scaffolding for which to build details later on. With that in mind, this primer isn’t a deep-dive tutorial or comprehensive guide. Instead, it aims to provide a high-level framework: what single-cell RNA-seq is, why it’s useful, the core steps of a typical workflow, and some common pitfalls to be aware of.

<hr>

### What Is Single-Cell RNA-Seq?

At its simplest, single-cell RNA-seq involves sequencing mRNA from individual cells by generating cDNA libraries tagged with barcodes that identify both the cell and each molecule sequenced. For the 10x Genomics Chromium platform, these barcodes are typically 16 bp (cell barcode) and 12 bp (unique molecular identifier, or UMI) in length.

This technology allows us to reconstruct a partial transcriptomic profile per cell, enabling insights into cellular identity, diversity, and state. While it may seem abstract, the basics of the wet-lab workflow help demystify how it works—and why certain analysis pitfalls can occur. Here’s a step-by-step overview (visual learners may find the field sketch below helpful as a guide):

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/Single_cell_diagram.png"
       style="max-width: 1400px;"
       alt="Single Cell Diagram" class="glow-frame">
  <figcaption>Wet-lab Single Cell Workflow</figcaption>
</figure>

##### Wet Lab Workflow

1. Tissue Resection
A physical tissue sample is obtained. For instance, we might collect a hepatocellular carcinoma specimen to study tumor-infiltrating lymphocytes (TILs) and the tumor microenvironment (TME).

2. Tissue Dissociation
The tissue is enzymatically and mechanically dissociated (often in a Falcon tube) into a suspension of viable, individual cells.

3. Microfluidics & GEM Formation
The cell suspension is loaded into a microfluidic chip, where Gel Bead-in-EMulsion (GEMs) are formed. Each GEM ideally contains one cell and one barcoded bead, encapsulated in oil. This creates thousands of tiny, isolated reaction chambers.

4. Library Preparation
Each bead is coated with oligos structured as follows:
[Gel-Bead - 5' - PCR Handle - Cell Barcode (16bp) - Random UMI (12bp) - Poly(dT) - mRNA 3']. 
Inside each GEM, cells are lysed, allowing mRNA molecules to bind to the poly(dT) tail. Reverse transcription produces single-stranded cDNA, droplets are broken, and PCR amplification follows. Finally, sequencing adapters (e.g., Illumina) are ligated to the DNA fragments.

5. Sequencing
Libraries are sequenced using standard sequencing platforms. The resulting data contain read sequences with embedded cell barcodes and UMIs.

<hr>

### Why Is This Useful?

Understanding what single-cell RNA-seq does is helpful—but the power lies in the questions it enables us to ask and answer. Here are two real-world examples from the field of cancer genomics to illustrate it's power:

1. Immune Checkpoint Inhibitor Response
Suppose you're studying a cohort of hepatocellular carcinoma patients treated with checkpoint inhibitors. Patients vary: some respond, others progress.
With scRNA-seq, you could ask:
- Are specific T-cell subtypes (e.g., exhausted vs. cytotoxic) enriched in responders/progressors?
- Can we identify biomarkers pre- or post-treatment which correspond to patient outcomes?
- Are tumor cells expressing markers of immune escape in progressive disease?

2. Tracking Tumor Evolution in an N-of-1 Case
In a patient receiving a novel therapy, early scRNA-seq timepoints might show response followed by resistance. By profiling transcriptomic changes over time, you could track sub-clonal evolution and transcriptomic response to treatment informing targeted changes to the treatment regimen.

While these examples focus on cancer, scRNA-seq is applicable across domains: microbiome studies, developmental biology, cardiac remodeling after infarction, and more. In fact, recent advances (as of 2020) now combine spatial data from histology with single-cell transcriptomics, adding positional context to cell expression profiles.

<hr>

### Core Analysis Steps

Once sequencing is complete, the computational pipeline begins. I focus here on the [10x Genomics platform](https://www.10xgenomics.com/platforms/chromium?utm_medium=search&utm_source=google&utm_term=10x+genomics+single+cell&useroffertype=website-page&utm_content=website-page&utm_campaign=701VI00000Fl8ZSYAZ&usercampaignid=701VI00000Fl8ZSYAZ&gad_source=1), however while specific tools may differ, nearly every scRNA-seq analysis follows the same conceptual steps:

1. Generate a Cell-by-Gene Count Matrix
- For 10x Genomics data, this is done using [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), which aligns reads and produces a matrix of counts:
(cells × genes) → count of transcripts per cell-gene pair.
  - In R, `Read10X()` from the Seurat package loads this matrix.
  - The matrix is usually sparse, meaning most values are zero. Tools like Seurat use efficient memory structures to omit explicit zeros.

2. Preliminary Quality Control (QC)
- Before analysis, remove problematic cells:
  - High mitochondrial gene content (>15–20%) may indicate dying or stressed cells.
  - Low or high total RNA counts may reflect empty droplets or doublets (more on that later).
  - Low gene complexity may suggest poor capture or degraded RNA.

3. Normalization
Different cells are sequenced to different depths. We must normalize to compare them fairly.
- A typical method:
Normalized Count = (UMI / Total UMIs per cell) × scale factor -> log-transformed
- This corrects for library size but assumes cells have roughly similar total RNA content. Alternatives like `SCTransform()` can address this more rigorously.

4. Feature Selection (Highly Variable Genes)
- We identify transcripts showing biologically meaningful variation, not just technical noise.
  - Low-expression genes with high relative variability may be stochastic.
  - Tools apply variance-stabilizing transformations (VST) to rank genes by informative variability.
*Optional: Rescaling can also be done to regress out confounders (e.g., batch effects).*

5. Dimensionality Reduction
- Reduce complexity while preserving key patterns:
  - Run PCA on variable genes.
  - Use an elbow plot to select significant principal components.

6. Clustering
- Group cells by similarity:
  - `FindNeighbors()` builds a graph based on PCA-reduced space.
  - `FindClusters()` identifies discrete communities (e.g., Louvain or Leiden algorithm).

7. Visualization (UMAP / t-SNE)
- Create a 2D layout (e.g., UMAP) of the high-dimensional data.
  - These visualizations preserve local structure and are useful for spotting patterns and clusters.

8. Annotation and Marker Discovery
- Now that we have clusters, we need to identify what they represent:
  - Run `FindMarkers()` (in Seurat) or equivalent to detect differentially expressed genes.
  - Use known markers (e.g., CD3E for T cells) to label clusters.

<hr>

### Common Pitfalls & Red Flags
Single-cell workflows are powerful, but not foolproof. Here are some things to watch out for:

1. Empty Droplets and Doublets
- Empty GEMs contain a bead but no cell -> no transcript capture.
- Doublets contain two cells -> one barcode, mixed transcriptome.
*Solution: Use QC filters (e.g., total UMI counts) and tools like DoubletFinder.*
*Refer to the figure below for a visual reference*

<div class="clearfix"></div>

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/emusion_beads.png"
       style="max-width: 1400px;"
       alt="Gel bead-in-EMulsion" class="glow-frame">
  <figcaption>Gel bead-in-EMulsion Conditions</figcaption>
</figure>

2. High Mitochondrial Content
- Elevated MT gene expression suggests the plasma membrane rupture or cell is undergoing apoptosis.
- Rates vary by tissue ([Mercer et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3160626/)), but >70% is a red flag.

3. Ambient RNA Contamination
- If MALAT1 or other housekeeping genes are ubiquitous across clusters, suspect ambient RNA.
- This arises from free-floating RNA (e.g., from lysed cells) being captured in the GEMs.
*Tool: SoupX can help estimate and subtract ambient signal.*

4. Flat or Indistinct UMAPs
- If your UMAP looks like a streak or blob, something went wrong.
- Common cause: over-normalization or failure to remove batch effects.

<hr> 

### Tools & Resources
- [Seurat](https://satijalab.org/seurat/) (R) – Most popular single-cell toolkit
- [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) (R) – More modular S4 framework
- [SoupX](https://github.com/constantAmateur/SoupX) (R) – Ambient RNA correction
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) - Detect doublet's in single-cell sequencing data
- [rnabio.org: Module 8](https://rnabio.org/module-08-scrna/0008/01/01/Intro_to_scRNA/) – Hands-on single-cell RNA-seq tutorial


