---
layout: post
title:  "scRNAseq - Primer"
date:   2025-10-09 10:13:50 -0500
description: "A Primer on Single-Cell RNA-seq: Concepts, Steps, and Pitfalls"
tags: [scRNAseq, bioinformatics]
---

### Introduction
Welcome to a primer on single-cell RNA sequencing (scRNA-seq). Studies suggest that new information is retained more effectively when paired with a high-level overview before diving into the details. With that in mind, this primer isn’t a deep-dive tutorial or comprehensive guide. Instead, it aims to provide a high-level framework: what single-cell RNA-seq is, why it’s useful, the core steps of a typical workflow, and some common pitfalls to be aware of.

### What Is Single-Cell RNA-Seq?

At its simplest, single-cell RNA-seq involves sequencing mRNA from individual cells by generating cDNA libraries tagged with barcodes that identify both the cell and each molecule. For the 10x Genomics Chromium platform, these barcodes are typically 16 bp (cell barcode) and 12 bp (unique molecular identifier, or UMI) in length.

This technology allows us to reconstruct a partial transcriptomic profile per cell, enabling insights into cellular identity, diversity, and state. While it may seem abstract, the basics of the wet-lab workflow help demystify how it works—and why certain analysis pitfalls can occur. Here’s a step-by-step overview (visual learners may find the sketch below helpful):

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/Single_cell_diagram.png"
       style="max-width: 1400px;"
       alt="Single Cell Diagram" class="glow-frame">
  <figcaption>Single Cell Workflow</figcaption>
</figure>

<hr>
##### Wet Lab Workflow

1. Tissue Resection
A physical tissue sample is obtained. For instance, we might collect a hepatocellular carcinoma specimen to study tumor-infiltrating lymphocytes (TILs) and the tumor microenvironment (TME).

2. Tissue Dissociation
The tissue is enzymatically and mechanically dissociated (often in a Falcon tube) into a suspension of viable, individual cells.

3. Microfluidics & GEM Formation
The cell suspension is loaded into a microfluidic chip, where Gel Beads-in-Emulsion (GEMs) are formed. Each GEM ideally contains one cell and one barcoded bead, encapsulated in oil. This creates thousands of tiny, isolated reaction chambers.

4. Library Preparation
Each bead is coated with oligos structured as:
[Gel-Bead - 5' - PCR Handle - Cell Barcode (16bp) - Random UMI (12bp) - Poly(dT) - mRNA 3']. 
Inside each GEM, cells are lysed, and mRNA binds to the poly(dT) tail. Reverse transcription produces single-stranded cDNA, droplets are broken, and PCR amplification follows. Finally, sequencing adapters (e.g., Illumina) are ligated to the DNA fragments.

5. Sequencing
Libraries are sequenced using standard Illumina platforms. The resulting data contain read sequences with embedded cell barcodes and UMIs.

### Why Is This Useful?

Understanding what single-cell RNA-seq does is helpful—but the power lies in the questions it enables us to ask and answer. Here are two real-world examples from the field of cancer genomics:

1. Immune Checkpoint Inhibitor Response
Suppose you're studying a cohort of hepatocellular carcinoma patients treated with checkpoint inhibitors. Patients vary: some respond, others progress.
With scRNA-seq, you could ask:
- Are specific T-cell subtypes (e.g., exhausted vs. cytotoxic) enriched in responders/progressors?
- Can we identify biomarkers pre- or post-treatment which correspond to patient outcomes?
- Are tumor cells expressing markers of immune escape in progressive disease?

2. Tracking Tumor Evolution in an N-of-1 Case
In a patient receiving a novel therapy, early scRNA-seq timepoints might show response followed by resistance. By profiling transcriptomic changes over time, you could track sub-clonal evolution and inform targeted changes to the treatment regimen.

While these examples focus on cancer, scRNA-seq is applicable across domains: microbiome studies, developmental biology, cardiac remodeling after infarction, and more. In fact, recent advances (as of 2024) now combine spatial data from histology with single-cell transcriptomics—adding positional context to cell epxression patterns.

### Core Analysis Steps

Once sequencing is complete, the computational pipeline begins. While specific tools may differ, nearly every scRNA-seq analysis follows the same conceptual steps:

1. Generate a Cell-by-Gene Count Matrix
For 10x Genomics data, this is done using Cell Ranger, which aligns reads and produces a matrix of counts:
(cells × genes) → count of transcripts per cell-gene pair.
- In R, Read10X() from the Seurat package loads this matrix.
- The matrix is usually sparse, meaning most values are zero. Tools like Seurat use efficient memory structures to omit explicit zeros.

2. Preliminary Quality Control (QC)
Before analysis, remove problematic cells:
- High mitochondrial gene content (>15–20%) may indicate dying or stressed cells.
- Low or high total RNA counts may reflect empty droplets or doublets (more on that later).
- Low gene complexity may suggest poor capture or degraded RNA.

3. Normalization
Different cells are sequenced to different depths. We must normalize to compare them fairly.
- A typical method:
Normalized Count = (UMI / Total UMIs per cell) × scale factor → log-transformed
- This corrects for library size but assumes cells have roughly similar total RNA content. Alternatives like SCTransform can address this more rigorously.

4. Feature Selection (Highly Variable Genes)
We identify transcripts showing biologically meaningful variation, not just technical noise.
- Low-expression genes with high relative variability may be stochastic.
- Tools apply variance-stabilizing transformations (VST) to rank genes by informative variability.
Optional: Rescaling can also be done to regress out confounders (e.g., batch effects).

5. Dimensionality Reduction
Reduce complexity while preserving key patterns:
- Run PCA on variable genes.
- Use an elbow plot to select significant principal components.
1. Running a principle component analysis on the data with features identified.
2. Identifying how many priniple components are potentially interesting, one method to do this would be an elbow plot examing the inflection point where we start to see diminishing returns.

6. Clustering
Group cells by similarity:
- FindNeighbors() builds a graph based on PCA-reduced space.
- FindClusters() identifies discrete communities (e.g., Louvain or Leiden algorithm).

7. Visualization (UMAP / t-SNE)
Create a 2D layout (e.g., UMAP) of the high-dimensional data.
- These visualizations preserve local structure and are useful for spotting patterns and clusters.

8. Annotation and Marker Discovery
Now that we have clusters, we need to identify what they represent:
- Run FindMarkers() (in Seurat) or equivalent to detect differentially expressed genes.
- Use known markers (e.g., CD3E for T cells) to label clusters.

### Common Pitfalls & Red Flags
Single-cell workflows are powerful, but not foolproof. Here are some things to watch out for:

##### Empty Droplets and Doublets
- Empty GEMs contain a bead but no cell → no transcript capture.
- Doublets contain two cells → one barcode, mixed transcriptome.
Solution: Use QC filters (e.g., total UMI counts) and tools like DoubletFinder.

<div class="clearfix"></div>

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/emusion_beads.png"
       style="max-width: 1400px;"
       alt="Gel bead-in-EMulsion" class="glow-frame">
  <figcaption>Gel bead-in-EMulsion</figcaption>
</figure>

#### High Mitochondrial Content
- Elevated MT gene expression (>20–30%) suggests membrane rupture or apoptosis.
- Rates vary by tissue (Mercer et al. 2011), but >70% is a red flag.

#### Ambient RNA Contamination
- If MALAT1 or other housekeeping genes are ubiquitous across clusters, suspect ambient RNA.
- This arises from free-floating RNA (e.g., from lysed cells) being captured.
Tool: SoupX can help estimate and subtract ambient signal.

#### Flat or Indistinct UMAPs
- If your UMAP looks like a streak or blob, something went wrong.
- Common cause: over-normalization or failure to remove batch effects.

6) After we have specific traits via our PCA, we can determine which cells are most similar to each other by clustering, this is accomplished by:
1. Build a network-graph based on trait similarity derived from the PCA. In essence determine how similar or disimillar cells are from one another forming groups.
2. Determine clusters, in other words once we know how similar or disimilar cells are define boundaries to categorize cells into distinct groups

7) At this stage you would typically create a UMAP or tSNE plot to visualize the now distnict groups of cells. If you are familiar with PCA, think of this as a bi-plot, but keep in mind that this is only two dimensions. It can be used to visualize the single -cell experiment however and is commonly used in publications.

8) As a final step we have identified clusters of cells corresponding to traits based off the transcriptomic profile, but we have yet to define them. We don't know what cluster of cells corresponds to T-cells or myeloid cells, what are healthy cells vs cancer. To do this we must identify specific biomarkers or hallmarks. What this means is performing almost a differential expression analysis, we are looking to see what transcripts, based on abundance estimates, define groups. We can then assign labels like "T-cells" to groups with inference or through mathematical means if we have a pre-defined list of biomarkers.

### Common pitfalls and Red flags

A complete breakdown of all pitfalls and ways to address QC issues is beyond the scope of this posts, but is perhaps a future topic. Instead I will mention a few nuances to be aware of. Let's start by going back to our standard single-cell workflow, specifically the microfluidics. In an ideal world a GEM would contain a single cell, with a single bead, this happens most of the time, but not always. It is possible for the emulsion bead, the combination of barcoded read and cell contained in it's own mini reaction chamber to diverage from this expectation. This is illustrated in the field sketch below, in a) we see an empty droplet, the barcoded bead is encapsulated but there is no cell, and so no usefull data will be generated. Conversley in c) we see another deviation from our expectation, the emulsion bead has encapsulated two cells, known as a doublet, this is problematic as both of these cells will share the same cell barcode, even though they could be fundamentally different. Additional combination, triplets for example, can form as well, fortunately these occurences are relatively rare and are easily removed as part of most pre-QC workflows. Specifically, total cell count cutoff's are generally applied, imagine a density plot of all reads for each cell, a doublet will have aprox. twice the number of reads as expected, and can be cut off using simple filters designed to capture that.

Mentioned in the preliminary QC step above, another pitfall relates to the percentage of MT DNA in a cell, while rates can vary amongst tissues, mercer et al. 2011 reports rates varying from 5-30%, rates should not be excessivley high. Greater than 70%, depending on the expected cells, would be cause for extreme suspicoun, and could be dead or dying cells, perhaps due to too much mechanical distruption while creating the cell suspensions. Whatever the cause, a high MT percentage indicates the cells plasma membrane has been compromised. Fortunately these can be easily removed by calculating the proportion of MT gene abundance in relation to the entire library and filtering out those cells which are outliers.

Of further concern would be the expression of MALAT1 or other "house-keeping" genes in every cluster of cells identified. This would be abnormal, and highly indicative of Ambient RNA contamination. This is a slightly different problem than what is previously described and indicates that free floating transcripts from dying cells were released into the cell suspension. These are then encapsulated in the emulsion bead and amplified along with the rest of the captured cells. I will briefly mention that there are tools that are designed to correct for this, one prominent one being the R library soupX, however optimization of the wet-lab protocol would be the ideal fix.

Lastly I will mention that distinct clusters should be observable in the UMAP plot in almost every conceivable instance. If there is a single cluster or one large streak, the problem is likley one of over-normalization and not that a single-cell type was sequenced.

### Tools & Resources
- Seurat (R) – Most popular single-cell toolkit
- SingleCellExperiment (R) – More modular S4 framework
- SoupX (R) – Ambient RNA correction
- rnabio.org: Module 8 – Hands-on single-cell RNA-seq tutorial


