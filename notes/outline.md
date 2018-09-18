
# AP2 gene expression regulates inflorescence architecture
> Working title

## Introduction

1. Why is inflorescence architecture / branching important (yield etc)? How is architecture different between domesticated and wild species?

2. How does gene expression regulate transitions in meristem identity? How do transitions in meristem identity relate to inflorescence architecture, yield?
    - "evolutionary diversity in Solanaceae inflorescence complexity is determined by subtle modifications of transcriptional programs during a critical transitional window of meristem maturation" [@Lemmon_Evolutioninflorescencediversity_2016]

3. How are Ap2 genes involved in meristem identity, architecture etc (maybe also stress)? (maybe in discussion)

> I agree to discuss AP2s (at least partially) in the introduction

4. Here we investigate the architecture of wild and domesticated rices panicles and their transcriptomes. We then focus and further characterize the AP2-EREBP family.

## Results

### 1. Methods summary

### 2. Panicle traits of domesticated Asian and African rice are similar. Panicles of domesticated species produce more primary and secondary branches and more spikelets.
- **Figure 1a**: photo of mature panicles (images panicle spread) Tom
- **Figure 1b**: Phenotyping results *e.g.* Cali phenotyping, MTP phenotyping, PCA (put asian domestic together) Tom
- correlation between spikelet number and branching complexity justifies transcriptome samples
- **SF 1** phenotype correlations.

### 3. There are "transcriptomic differences" between [domesticated and wild] varieties.
- **Figure 2**: global view of transcriptome *e.g.* PCA (Otho)
- **SF 2**: Montpellier phenotyping (Tom)
    + re-order phenotypes
- **SF 2b**: panicle samples
- **SF 3**: distance heatmap (Tom)
- **SI**: DE genes table
- Topology / mathematics
- african species cluster together, asian species cluster together -> domestication changes are minor compared to species-level changes

> Comment OM: I wouldn't stress "yield" too much

### 4. AP2-EREBP are the most prominently DE family in transcriptomes of developing panicles.
- **Table 1**: TFs differentially expressed between stages in all species (featuring AP2 genes). Try z-score for AP2s (Tom).
- **Figure 3a**: PCA enrichment + AP2 amd HB heatmap (Otho).
- PCA splits genes expressed in BM and SM on component 5, AP2/ERF are the most enriched across PC5
- at least 20 AP2 are among the 400 genes preferentially expressed in BM, at least.
- DREB ERF and RAV peak earlier, and AP2 peaks later
- Can mention HB / NAC, other families here
- **SF 4** PCA enrichement + MADS, NAC, MYB.
- **SF 5** PCA enrichement all families.
- [Differential] Expression of AP2 genes can be grouped by subfamily **Figure Supp  5**: AP2-EREBP phylogeny with gene expression (Otho + My) - **It would be nice to have a phylogeny of all ap2s**


### 5. AP2-EREBP expression is associated with high yielding indica.
- **Figure 4**: clusters and correlation with phenotypes (or PCA?) (Tom)
- **Figure 5b / Table 3**. DE genes for interaction between stage and species (African and Asian) (Tom). (enriched for AP2s).
- **SF**: Cluster 5 only
- AP2s in cluster 4 & 5, generally they behave opposite in indica, some of them also in PC5
- Cluster1: Arag1; Cluster2: EREBP86; Cluster4: PLT9, IDS1, EREBP153; Cluster5; AP37/ERF3, DREB4-1, DREB4-2

### 6. AP2-EREBP genes are differentially expressed over a timecourse of panicle development:
- **Figure 5**: fluidigm results for AP2 (Otho)
- If this only confirms 4., we could move it to supp. If the different stages add something new, keep it in the main paper.

### 7. Differences in AP2 gene (promoter regions+CDS) are associated with differences in expression
- **Figure 6**: (Otho and Hélène will try this experiment)

### 8. AP2 mutants have defects in panicle branching.
- **Figure 7**: mutant characterisation (Otho / Hélène)
- Also check LMD data for these (Tom). There's already an AP2 figure from the LMD data, Tom will find it and send it out.
- PLT8 TPM in supplementary

**Mutant summary:**

**LOC_Os08g31580 / erf48**

- Typical expression pattern of PC5, high ranking, mostly expressed in BM, (Generally, does indica starts developing earlier?)
- DREB gene,
- effect on RACHIS and PRIMARY branches:
    - Rachis shorter,
    - More Primary Branches (Maybe???),
    - Primary branches shorter,
    - Primary branch internodes shorter,
    - No effect on secondary branches,
- Also less variance in all features in the mutants

**LOC_Os07g03250 plt8 / crl 5**

- DE only in indica, Higher in BM,
- ERF gene,
- Always less variance in the mutant,
- Effect on Rachis and primary branches:
    - Rachis shorter,
    - Less primary branches,
    - Both order branches longer?
    - Maybe more secondary branches and more spikelet? Small effect,

**LOC_Os05g32270 erf142 &&& LOC_Os06g03710 smo2 ???**

They behave always the same, erf142 has more extreme values.

- erf142 classified as DREB by Sharoni et al and as Unusual AP2 by Aya et al, 2014, which is the one that reported the mutant.
- Everything is smaller, besides, strangely, SbIntL, which gets longer.
-  It makes significantly less Sb and Sp, which fits with its function in Auxin signalling pathway.



## Discussion

1. We found lots of evidence that AP2 genes regulate panicle branching. Discuss different mechanisms *e.g.* repression of determinate state, promotion of indeterminate state, control exhaustion of branch meristem and identity / rate of production of axillary meristem and identity of spikelet meristem. Links to stress response. Discuss potential for changes in AP2 gene expression to be involved in evolution of branching phenotype.

2. Limitation of these results: lack of functional analysis, mutants etc. Potential for redundancy between genes. Potential contributions by genes from other families to regulatory modules.

3. How could this apply to [crop] domestication in general? What are the similarities in other branching processes e.g. root development, inflorescence development in other species?

4. What does the AP2 story us about the evolution/development of branching architecture?

## Notes

- What about the ARF genes?

> Note OM: While AP2s stand out, we have nice observations on many other families. Shall we draft the paper fully on AP2 and then try to add in the most interesting observations about the other families and see how they fit?
