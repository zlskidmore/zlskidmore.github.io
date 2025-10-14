---
layout: post
title:  "scRNAseq - Primer"
date:   2025-10-09 10:13:50 -0500
description: "A Primer on Single-Cell RNA-seq: Concepts, Steps, and Pitfalls"
tags: [scRNAseq, bioinformatics]
---

### Introduction

Welcome to a primer on single-cell RNAseq. I've read that information is retained more effeciently if a high level summary or general knowledge is available prior to diving into a topic. As such this primer is not intended to explain every nuance on single-cell analysis nor is it intended to be a tutorial, rather I will focus on a high level what single-cell sequencing is, why it is usefull, the core steps in all single-cell analysis studies, and common pitfalls. 

### What is single-cell sequencing?

So what is single-cell RNAseq? at the simplest level it is sequencing of mRNA, from cDNA libraries with the addition of sequence barcodes in order to identifiy both individual cells and individual molecules, for the 10x genomics kits these barcodes are 16bp and 12bp in length respectively. These barcodes are incorporated into the sequenced transcripts allowing for us to map a partial transcripomic profile to individual cells. It's amazingly impressive technology, but first, it is usefull to understand the wet lab workflow, at least on a surface level to understand not only how the core technology works but some pitfalls that we can get into later. To that end the basic steps are as follows, if you are a visual learner like me please refer as well to the journal sketch below:

1. Tissue Resection: The tissue of interest is resected to get an entire piece of tissue of interest. This could for example be a hepatocellular carcinoma and we are interested in the makeup of cells from the tissue, particularly any TIL's.

2. Tissue Dissociation: The solid tissue of interest is broken down into a suspension of individual viable cells, generally done in a falcon tube it involves enzymatic digestion and mechanical distruption in order to isolate cells.

3. Microfluidics: The Tissue suspension in step 2 is run through microfluidic channels. This is the core of single-cell technology and essentially involves encapsulating single cells with barcoded beads using oil, which forms a sort of individual reaction chamber, known as a single-cell GEM. That is to say each GEM formed in oil contains in theory one cell, and one barcoded gel bead.

4. Library Prep: Each bead is covered in oligos, and take the following form: [Gel-Bead - 5' - PCR Handle - Cell Barcode (16bp) - Random UMI (12bp) - Poly(dT) - mRNA 3']. Within each droplet cells are lysed and the mRNA covalently bonds to the Poly(dT) tail, which is already attached to the bead. Reverse transcriptase enzyme is introduced which forms single stranded cDNA and droplets are subsequently broken. PCR is then performed to amplify cDNA fragments producing double-stranded cDNA. Next adapters (i.e. illumina) are ligated on to allow for fragments to bind to the flow-cells on a sequencer.

5. Sequencing: Resulting libraries can then be sequenced on a sequencer per the ussual protocols which are out of scope of this single-cell RNAseq primer.

### Why is this usefull?

You should now understand in principle what single-cell does, but why should we care? what does this allow? To answer that i'm going to give a couple of hypothetical situations from my own background in cancer genomics:

1) Let's imagine you are sequencing a cohort of hepatocellular carcinoma patients treated with an immune checkpoint inhibitor. We have observed clear outcomes from response to therapy to progressive disease and we want to know what in the TME is differentiating the responders v non-responders.
a) are there spcific T-cell subpopulations between the groups for which we can differentiate cytotoxic or exhausted T-cells? If so how many cells are in each population?
b) assuming we had single-cell data before/after treatment, are there specific biomarkers which can be differentiated in responding patients, perhaps a sub-population of tumor cells have found a mechanism for immune escape.

2) As an additional example, perhaps we have a N-of-1 case study being treated with a novel therapy. Over time we see an initial response followed by progression for which we have single-cell data for. We can perhaps track the transcriptomic profiles of known tumor sub-clones in order to gain a deeper understanding of this hypothetical patients sub-clonal tumor interactions with the therapy. Perhaps with some luck we can even make an informed modification to the treatment regimine to target multiple expanding sub-clones.

These are just some simple examples relevant to cancer genomics, however one could also imagine applications regarding microbiome populations as it relates to gut health or understanding individual cellular interactions after a myocardial infarction. In fact the technology, relatively recently as of 2024, has been pushed even further incorporating spatial data (i.e. from histological slides) in conjuction with single-cell to gain an even deeper understanding.

### Core steps

So after we have our sequencing data, what are the core steps analysis steps. They can be cleanly divided, in my opinion, into 5 distinct steps. I outline these steps below, and while this primer is not intended to explicitly state how to do each of these steps, they should be recognizable in almost every single-cell workflow, regardless of the underlying analysis libraries used.

1) Convert the raw sequencing data to a count based matrix where the matrix is unique cell barcodes x Transcripts with cell values representing cell counts. For the 10x Chromium platform there is a specific pipeline known as Cell Ranger which will achieve this. After the transcriptome count matrix is generated it can be loaded and analyzed, depending on the tools used there may be specific convenience functions for this, for example the R library "Seurat" has a function called `Read10x()` specifically for output from 10x genomics Cell Ranger.

a) As an interesting aside, you will notice that typically the matrix is sparse, in other words for a typical experiment only a fraction of the transcritome is captured for each cell resulting in many 0 counts. Instead of encoding 0's which takes up space in memory, nothing is encoded for those values and tools like Seurat will print out place-holders for those matrix cells instead. A bit more computationally intensive, but well worth it when taking into account the memory savings.

2) Preliminary QC is performed. This can become more advanced but essentially we want to remove those cells which appear on the surface to have problems. We actually go into this in a bit more detail in the next section however consider this as some of the cells/data we might wish to remove:

1. Dying cells or otherwise non-viable cells, these can typically be identified by high concentrations of MT DNA
2. Cells with abnormal counts, these outliers could represent problems with the reaction chamber which is described above, known as GEMs.

3) Matrix count normalization is performed, in other words each captured and sequenced cell needs to become comparable to the other captured cells. If you are familiar with traditional RNAseq experiments think of this as controlling for library size. We need to make sure that a cell doesn't have higher abundance counts for a gene just because it was sequenced deeper. A simple approach commonly used is a log-normalization, in other words something like (cell gene reads / sum of cell gene reads) × scale factor -> log transform. This however does assume each cell has aprox. the same number of molecules. There are ways around this however if that assumption is violated.

4. Identify initial features, i.e. transcripts of interest. What I mean by this is that we need to pull out suspected biological signal vs noise. As an example a transcript can have a high variation amongst cells, however if the counts for that transcript in each cell are low this could not be biology but rather simple noise due to lower general transcript abundance. Put another way consider a simple but extreme example, your sequence library does not contain millions of reads but only 100, and two of those reads happen to map to the same gene, for one cell but nothing else for any other cell. There are bigger problems if this ever truely happens of course, but how do we know if those two reads are significant, that is what we are trying to accomplish here. A common method to do this is known as a variance stabilizing transform or VST.

4a). This step is somewhat optional, but commonly done so is included here, after a VST or other method of identifying features of interest we re-scale the data. In other words we are taking out or at least limiting the effect other covariates may have such as batch.

5) Perform a dimensionality reduction. What we want is to take our now cleaned and normalized count data and determine for that data what traits define groups. Recall that while we have data for each cell, we don't which traits define each group of like cells or know which cells are most similar to each other. A dimensionality reduction will clear the first part of that question. This is typically done by:

1. Running a principle component analysis on the data with features identified.
2. Identifying how many priniple components are potentially interesting, one method to do this would be an elbow plot examing the inflection point where we start to see diminishing returns.

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

### Tools and additional resources

- Seurat, single-cell analyis library in R, excellent vignettes and workflows are available
- SingleCellExperiment, A Bioconductor package, a bit more flexible than seurat at the cost of added complexity
- soupX, R library designed to remove cell free mRNAs from single experiments
- rnabio.org - Module 8, a hands on tutorial for a single cell experiment






Lorem ipsum dolor sit amet consectetur adipiscing elit. Quisque faucibus ex sapien vitae pellentesque sem placerat. In id cursus mi pretium tellus duis convallis. Tempus leo eu aenean sed diam urna tempor. Pulvinar vivamus fringilla lacus nec metus bibendum egestas. Iaculis massa nisl malesuada lacinia integer nunc posuere. Ut hendrerit semper vel class aptent taciti sociosqu. Ad litora torquent per conubia nostra inceptos himenaeos.
<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/Single_cell_diagram.png"
       style="max-width: 1400px;"
       alt="Single Cell Diagram" class="glow-frame">
  <figcaption>Single Cell Workflow</figcaption>
</figure>
Lorem ipsum dolor sit amet consectetur adipiscing elit. Quisque faucibus ex sapien vitae pellentesque sem placerat. In id cursus mi pretium tellus duis convallis. Tempus leo eu aenean sed diam urna tempor. Pulvinar vivamus fringilla lacus nec metus bibendum egestas. Iaculis massa nisl malesuada lacinia integer nunc posuere. Ut hendrerit semper vel class aptent taciti sociosqu. Ad litora torquent per conubia nostra inceptos himenaeos.


<div class="clearfix"></div>

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/scRNAseq_primer/emusion_beads.png"
       style="max-width: 1400px;"
       alt="Gel bead-in-EMulsion" class="glow-frame">
  <figcaption>Gel bead-in-EMulsion</figcaption>
</figure>

