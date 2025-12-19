---
layout: post
title:  "Mimimal Residual Disease In Cancer"
date:   2025-12-10 10:13:50 -0500
description: "MRD in cancer, where we are and current trends"
tags: [MRD, liquid_biopsy]
---

# Minimal Residual Disease: What It Is, and How We Approach It

<hr>

### What is minimal residual disease?

Put simply, minimal residual disease (MRD) refers to cancer that appears to be gone after treatment, yet persists at levels too low to detect with standard clinical tools. Imaging may show a complete response, blood counts may normalize, but a small population of malignant cells can remain and eventually spark a relapse.

To ground this in a real example, consider a patient with small cell lung cancer (SCLC). These tumors almost universally carry TP53 and RB1 mutations and are typically treated with a platinum-based chemotherapy. Platinum agents induce double-strand breaks, and because tumor cells divide rapidly, they generally take the brunt of the damage. Many patients show an initial response.

But curing cancer requires eliminating every malignant cell — and that’s the hard part.

Tumors are rarely uniform. They evolve rapidly, generating genetically distinct subclones with different levels of drug sensitivity. This phenomenon, tumor clonality, is one of the central reasons early detection matters. The more heterogeneous a tumor becomes, the more likely it is that at least one subclone can survive therapy.

So even when a cancer appears “cleared” on imaging, it may not be gone. A few surviving cells can continue dividing quietly until clinical relapse becomes inevitable. In SCLC this pattern is common: patients often look like they’ve responded well, yet molecularly the disease persists, and recurrence can happen with startling speed.

That lingering, invisible population of cells — too small for a scan to detect but biologically active — is MRD.

<hr>

### Why are we developing MRD assays?

With a working definition of MRD, we arrive at the next question, why do we care? The short answer is that if we can detect cancer while it still exists only at a molecular level, we can intervene earlier and more effectively. Knowing that disease is still present allows us to either extend therapy, escalate treatment, or switch treatment modalities before the tumor has the opportunity to rebound and diversify.

A central truth in oncology is that time favors the tumor. The longer cancer cells are allowed to persist, the more opportunity they have to evolve into treatment-resistant subclones. Catching a molecular relapse early, weeks to months before a CT scan can show anything, gives a meaningful head start to head off further clonal evolution.

Beyond imrpoved treatment outcomes, a practical advantage also exists. Blood-based MRD assays are far more accessible than imaging. Not every patient lives near a center with radiology capabilities and a ct-scanner, delays in scheduling or geography can push imaging out by critical weeks. A blood draw is simple. It’s cheap.

That combination, earlier detection, convenience, and wider accessibility, is why the field is investing so heavily in MRD technologies.

<hr>

### Why is ctDNA an attractive target for MRD assays?

We’ve discussed what MRD is and why it matters, so the next question is, how do we actually detect it? Most modern MRD assays rely on a liquid biopsy, a blood draw is collected in a Streck or EDTA tube and then analyzed molecularly. Other sample types exist in theory, and I’ll touch on them later, but ctDNA remains the backbone of almost every MRD assay in development today.

The reason is sensitivity. MRD requires us to detect cancer at astonishingly low levels of signal. We’re often looking for tumor-derived DNA (ctDNA) fragments that represent 0.1% of the total cell-free DNA, even an order of magnitude below that is preferable. Achieving that level of detection isn’t just a technical challenge; it narrows the viability of the biological materials we can work with.

ctDNA hits core checkmarks offering a rare combination of attributes:

- It’s accessible: a blood draw is easy to obtain compared to specialized imaging
- It’s cost-effective: a blood draw can be processed and sequenced relatively quickly and easily
- It’s biologically informative: ctDNA carries tumor-specific biomarkers that are detectable and measurable
- It’s scalable: NGS sequencing pipelines can process thousands of samples in a lab each week

All of this while also approaching useful limit of detection and limit of blank thresholds to matter.

<hr>

### Key MRD approaches, strengths and weaknesses

##### Tissue Informed Methodolgies

Many of the leading MRD assays are tissue-informed, which pushes the field squarely into the realm of personalized medicine. The basic idea is fairly straightforward, you sequence a patient’s tumor tissue up front, identify the somatic mutations that define it, particularly the founding variants, and then use that information to inform where you look in the liquid biopsy.

Once the tumor-specific variants are known, you can build a custom capture panel targeting specifically those sites and sequence the patient’s plasma to an extreme depth, often in the range of 1,000–10,000×. These assays typically incorporate unique molecular identifiers (UMIs), which tag each original DNA fragment before any PCR amplification occurs. When multiple reads originate from the same tagged fragment, you can then collapse them into a consensus sequence, dramatically reducing the chance that a random PCR error is classified as a true somatic variant.

On a first pass this might seem like overkill, however MRD detection is fundamentally a needle-in-a-haystack problem. We’re looking for the faintest of signals in a background dominated by healthy cell-free DNA. Anything that reduces noise or boosts signal confidence helps.

The strengths of tissue informed approaches are fairly intuitive:
- High sensitivity: You’re only looking for variants you know are present in the tumor
- High specificity: CHIP variants (more on these later) and other sources of background noise are largely filtered out
- Structured targeting: Ultra-Deep sequencing is concentrated at the most informative loci rather than spread across the genome

But as with any approach there are trade-offs. Critically, tissue informed MRD requires access to tumor tissue, something which might not always be available, especially in cancer types where biopsies are risky or yield very little material. Even when tumor tissue is available, sequencing, analysis, and panel design add time, cost and complexity. The extra logistical steps slow assay turnaround and make the approach less flexible than tissue agnostic alternatives.

##### Tumor Agnostic Methodolgies

While there is no question, tissue-informed approaches currently offer the best sensitivity, the requirement for tumor tissue is a real limitation and drawback due to the reasons discussed above. Because of this, several creative and innovative MRD strategies have emerged that rely purely on plasma, without any aprior knowledge of a tumors mutational profile.

Let’s walk through a few of these.

###### Variant Based MRD

The absence of a personalized tumor panel doesn’t mean variants are off the list of viable biomarkers, it just makes the problem more difficult. Many cancers carry recurrent hotspot mutations. As I mentioned earlier, SCLC often carries TP53 and RB1 mutations; colorectal cancers (CRC) commonly harbor mutations in KRAS, NRAS, or BRAF. In a tumor-agnostic approach, we can look broadly for mutations like these rather than targeting a defined set of patient-specific personalized variants.

The challenge is that without ground truth, a gold standard, we’re more prone to misclassifying technical artifacts as somatic mutations, especially given the extremely low allele fractions involved in MRD. Lets consider a CRC case in which we detect a low-level KRAS mutation. Most sequencing workflows rely on PCR during library prep, and polymerases occasionally introduce errors, a mistaken base substitution that can look like a real somatic mutation. UMIs help here by tagging each original fragment before amplification. Unless an error happens in the very first PCR cycle, which can occur, UMI consensus logic will collapse the amplified reads back to a correct consensus sequence, dramatically reducing this type of issue.

Sequencing errors from the machine are another source of noise. Despite the very high accuracy, a sequencer processing millions of DNA fragments, and will occasionally misread a base, calling an C a T for example. Duplex sequencing offers a solution to control for this, because true mutations exist on both the forward and reverse strands of the DNA, a mutation that appears on only one strand of a DNA fragment is suspect. Duplex isn’t cheap, but it is a clean way to suppress these types of technical artifacts.

We must also consider the biology itself. We would be remiss not to mention CHIP variants (clonal hematopoiesis of indeterminate potential). As people age, hematopoietic stem cells will accumulate mutations and expand into small clonal populations. These mutations appear in plasma even though they have nothing to do with the tumor, often at similar levels to MRD detection thresholds. Fortunately these can be controlled for by sequencing the buffy coat (essentially a WBC control) alongside the plasma. In this manner most CHIP variants can be identified and excluded. This does come at the expense of a separate library prep and extra cost.

Finally I’ll end this section with a quick note on phased variants. Across all these approaches the goal remains the same, either reduce noise or boost signal, ideally both. Phased variants, which are two or more mutations occurring close enough together to be captured on the same DNA fragment, are a particularly interesting way to do this. The probability of two independent sequencing or PCR errors occurring in the exact same molecule are.... shall we say improbable, so identifying phased variants provide an extra layer of confidence when calling a variant as tumor derived.

###### Copy Number and LOH MRD

In the context of MRD, another source of signal beyond simply point mutations is aneuploidy, broad gains, losses, and structural variations in the cancer genome. Tumors can be remarkably chaotic. One region may show the deletion of an entire chromosome arm, while another contains high-level focal amplifications. Some may even exhibit copy-neutral loss of heterozygosity (LOH) where a region is deleted on one chromosome, but the remaining allele is then duplicated, returning the locus to a diploid copy number while still losing heterozygosity.

These events are, in some ways, easier to detect than single nucleotide variants. Instead of relying on a change at a single base, the signal is distributed across a much larger region of the genome. That broader footprint means fewer sequencing reads are needed to see the signal emerge, and the noise from individual basepair errors is far less relevant.

The trade-off here is sensitivity. Because the signal is spread out, aneuploidy-based MRD generally cannot reach the same detection limits as variant-based approaches. A practical example is ichorCNA, a widely used tool that estimates ctDNA fraction using low-coverage whole-genome sequencing. It works very well for moderate tumor fractions, but given the need for whole genome sequencing (WGS), and as such a restricted sequencing depth, the 0.1% thresholds we care about for MRD are elusive.

As with any MRD modality, the signal-to-noise ratio is everything. Copy-number inference brings its own technical artifacts, GC bias, mappability differences, etc. all need to be addressed before the biological signal becomes interpretable.
The broader point to make is that tumor-derived signal exists in many forms. Somatic variants are powerful, but they aren’t the only lens through which MRD can be viable.

###### Fragmentome MRD

Fragmentomics is the third major feature set used in liquid biopsy based MRD, and worth some attention. One of the foundational observations is that fragment sizes differ between normal cell-free DNA and tumor-derived cfDNA. Tumor fragments tend to skew shorter, and this size shift creates a detectable biomarker in plasma samples.

Fragmentomics extends beyond fragment length however. Tumor derived DNA often carries distinct end-motif patterns, reflecting altered nuclease activity in cancer cells, DNASE1 being a common example. These biases create reproducible signatures that can help distinguish ctDNA from normal cfDNA.

A further layer of interest comes from nucleosome positioning. If you look at coverage around transcription start sites (TSS), transcription factor binding sites (TFBS), or the first exon–intron junctions, you will often see characteristic dips in read depth. The logic is straightforward, DNA wrapped tightly around nucleosomes is physically protected from degradation in the bloodstream, while DNA that becomes unwound, for example, during active transcription, is more vulnerable and breaks down more readily. Tumor cells with abnormal chromatin states produce can distinctive patterns of nucleosome depletion, which show up as these dips in coverage at specific sites.

This creates a kind of chromatin fingerprint, it's not quite gene expression in the traditional RNAseq sense, but an adjacent signal that tumors leave will leave behind.

As always, there are trade-offs. Fragmentomics typically requires whole-genome sequencing, even if at shallow depth, which shifts the balance between sensitivity and cost. You gain broad structural information but sacrifice the ultra-deep coverage achievable with targeted variant assays.

Still, this is another orthogonal signal class, fragmentomics has real value. It captures tumor biology beyond that of variants alone cannot, and as such is worth mentioning in our MRD discussion.

<hr>

### Upcoming or Alternative Technologies

I want to briefly mention some alternative technologies and avenues of exploration. MRD detection isn’t limited to ctDNA, even though that’s where most of the focus sits today. As an example, cancer cells shed extracellular vesicles (EVs) which contain tumor-derived DNA, RNA, and proteins. These vesicles can be isolated with microfluidic systems and, in theory, provide a source of tumor specific biomarkers. The challenge here is scalability, current EV platforms are elegant but currently lack the throughput needed for widespread adoption.

Another direction I find particularly interesting is the idea of looking not for the tumor itself, but for tumor-adjacent signal. I am speaking of course about the immune system, even if the immune system can’t eradicate the cancer, it should still see and responds to it. If we reach a point where we can reliably map T-cell receptor (TCR) clonotypes to specific neoantigens presented on patient MHC molecules, then shifts in the TCR repertoire could, in principle, serve as an indirect MRD readout. Instead of tracking the tumor’s signal directly, we could track the immune system’s memory of it.

<hr>

### Final Thoughts

At the end of the day, MRD detection boils down to two questions:

- How do we boost signal?
- How do we reduce noise?

On paper, that sounds almost trivial. In practice, as I hope this overview has made clear, it’s far more complex. Every modality:

- variants
- copy number
- fragmentomics
- other experimental avenues

runs into its own set of biological and technical limitations.

My view is that multi-omic MRD assays are the next logical step. No single signal class is perfect, but each brings something complementary. Variants offer specificity, aneuploidy and structural features add robustness, fragmentomics broadens sensitivity, and emerging technologies may fill in remaining gaps. When combined, these layers can push detection limits low enough to make MRD monitoring mainstream and genuinely transformative for patient care. Indeed I think we may be on the cusp of this already.
