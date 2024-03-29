
# AP2 & HB gene expression pattern associated with inflorescence architecture diversity

> Working title (too vague IMO)
 
## Introduction

1. Why is inflorescence architecture / branching important (yield etc)? How is architecture different between domesticated and wild species?

2. How does gene expression regulate transitions in meristem identity? How do transitions in meristem identity relate to inflorescence architecture, yield?
    - "evolutionary diversity in Solanaceae inflorescence complexity is determined by subtle modifications of transcriptional programs during a critical transitional window of meristem maturation" [@Lemmon_Evolutioninflorescencediversity_2016]

3. How are Ap2 genes involved in meristem identity, architecture etc (maybe also stress)? (maybe in discussion)

> I agree to discuss AP2s (at least partially) in the introduction

4. Here we investigate the architecture of wild and domesticated rices panicles and their transcriptomes. We then focus and further characterize the AP2-EREBP family.

## Results

### 1. Methods overview

### 2. Panicle traits of domesticated Asian and African rice are similar. Panicles of domesticated species produce more primary and secondary branches and more spikelets.
- **Figure 1A–D**: photo of mature panicles (images panicle spread) Tom
- **Figure 1E–F**: Phenotyping results *e.g.* Cali phenotyping, MTP phenotyping, PCA (put asian domestic together) : just present PC1. Note the sequenced accessions.
- **SF 1** : other components

### 3. There are "transcriptomic differences" between [domesticated and wild] varieties.
- **Figure 2**: global view of transcriptome *e.g.* PCA (Otho).
    - Improve colors
- **SF 2**: Montpellier phenotyping (Tom)
- **SF 3** : phenotypes from Mtp samples (correlation between spikelet number and branching complexity justifies transcriptome samples)
- **SF 4**: panicle samples + Samples validation by QPCR using LHS1 & TAW
- **SF 5**: distance heatmap (Tom)
- **Table S1**: Exonic mapping / qc
- **Table S2**: all DE genes
- Topology / mathematics : Figure prepared by helene

- african species cluster together, asian species cluster together -> domestication changes are minor compared to species-level changes
- description of the dataset of each PC..and the finish by PC5

> Comment OM: I wouldn't stress "yield" too much

### 4. AP2-EREBP are the most prominently DE family in transcriptomes of developing panicles.

Global to TF enrichment : supp data

- **Table 1**: TFs differentially expressed between stages in all species (featuring AP2 genes). Try z-score for AP2s (Tom).
- **Figure 3**: PCA enrichment + AP2 amd HB heatmap (Otho) : justification of HB presentation and nopt MADS/NAC or other most statistically strong
- PCA splits genes expressed in BM and SM on component 5, AP2/ERF are the most enriched across PC5
- at least 20 AP2 are among the 400 genes preferentially expressed in BM, at least.
- DREB ERF and RAV peak earlier, and AP2 peaks later
- Can mention HB / NAC, other families here
- **SF 4** PCA enrichement + MADS (add subfamillies), NAC, MYB.
- **SF 5** PCA enrichement all families.
- [Differential] Expression of AP2 genes can be grouped by subfamily **Figure Supp  5**: AP2-EREBP phylogeny with gene expression (Otho + My) - **It would be nice to have a phylogeny of all ap2s (phylogeny of subclass in each species will be done later)**

### 5. AP2-EREBP expression is associated with high yielding indica.
- **Figure 4**: clusters and correlation with phenotypes (or PCA?) (Tom) / on the figure, remove r-enrichment fig and place just the name
- **Figure 5b / Table 3**. DE genes for interaction between stage and species (African and Asian) (Tom). (enriched for AP2s).
- **SF 6**: Cluster 5 only
- AP2s in cluster 4 & 5, generally they behave opposite in indica, some of them also in PC5
- Cluster1: Arag1; Cluster2: EREBP86; Cluster4: PLT9, IDS1, EREBP153; Cluster5; AP37/ERF3, DREB4-1, DREB4-2

### 6. AP2-EREBP genes are differentially expressed over a timecourse of panicle development:
- **Figure 6**: fluidigm results for AP2 (Otho) : remove japonica
- If this only confirms 4., we could move it to supp. If the different stages add something new, keep it in the main paper.
- What are the conclusions?

### 7. Differences in AP2 gene (promoter regions+CDS) are associated with differences in expression
- **Figure 7**: (Otho and Hélène will try this experiment) : **RAV2** duplication/check in 3000genomes.
  - More details: Can we detect which read comes from the each of the two duplicated genes? Can we detect this duplication in other indicas?
- **LAX1** also?

### 8. AP2 mutants have defects in panicle branching.
- **Figure 8**: mutant characterisation (Otho / Hélène)
- Also check LMD data for these (Tom). There's already an AP2 figure from the LMD data, Tom will find it and send it out.
- PLT8 TPM in supplementary.
- Remove erf48. Show gene expression in fluidigm and rnaseq.

## Discussion

1. We found lots of evidence that AP2 genes regulate panicle branching. Discuss different mechanisms *e.g.* repression of determinate state, promotion of indeterminate state, control exhaustion of branch meristem and identity / rate of production of axillary meristem and identity of spikelet meristem. Links to stress response. Discuss potential for changes in AP2 gene expression to be involved in evolution of branching phenotype.

2. Limitation of these results. Limitations of mapping against japonica, and strategies we used to deal with them, e.g. only looked at L2FCs not raw counts/expression, PCA, etc. Lack of functional analysis, mutants etc. Potential for redundancy between genes. Potential contributions by genes from other families to regulatory modules.

3. How could this apply to [crop] domestication in general? What are the similarities in other branching processes e.g. root development, inflorescence development in other species?

4. What does the AP2 story us about the evolution/development of branching architecture?

> some graphical model of the [parallel] changes in panicle architecture during domestication?

## Notes

- What about the ARF genes?

> Note OM: While AP2s stand out, we have nice observations on many other families. Shall we draft the paper fully on AP2 and then try to add in the most interesting observations about the other families and see how they fit?
