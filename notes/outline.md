
## Introduction

* evo devo and gene expression
* transcription factors / evolution of gene regulation
* domestication and gene expression (maize, tomato transcriptomes)
* details about rice domestication
* **aims**:
    - changes in phenotype, gene expression caused by / associated with domestication
    - common genes that respond to selection in different domestications
    - evolution of transcription factor family expression
    - does selection for the same phenotype affect expression of the same genes / sets of genes

## Results/discussion

**hypothesis**

- there is a group of genes that are repressed at the SM stage in wild species (-ve L2FC PBM vs. SM). These genes "promote indeterminate state". These genes are *not* repressed in the domesticated species, therefore the indeterminate state is maintained for longer.

**evidence**

1. Phenotype-L2FC correlation

- not sure yet...

2. PC5

- This component splits the two stages and it is mostly driven by the wild species: while nipponbare and especially indica are much closer to 0.
- The wild species display a much stronger burst of MADS-box related to differientiation
- The AP2 story might be more complicated, but indeed many of them are strongly down-regulated in most species, but they behave differently in indica, nipponbare

3. clustering of TF genes

- note that we are not necessarily limited to family-level comparisons. e.g. AP2 could regulate MADS etc.
- Hélène is currently looking for biologically significant stuff in the clusters
- no family survives bonferroni correction in cluster 6
- AP2 genes are enriched in cluster with high +ve L2FC in indica

4. differential expression

- note non-overlap of glab / indica genes. Of course selection could have acted on different genes in different domestications
- genes with a +ve L2FC for the stage:species interaction fit this category

### 1. Phenotypic analyses of African and Asian domesticated and wild rice

* Parallel evolution of phenotype during African and Asian domestications
* Correlations between *e.g.* spikelet number and branching complexity
    - these correlations justify the choice of PBM and SM for sequencing
    - also, "evolutionary diversity in Solanaceae inflorescence complexity is determined by subtle modifications of transcriptional programs during a critical transitional window of meristem maturation" [@Lemmon_Evolutioninflorescencediversity_2016]
* **Figure 1**:
    - phenotyping 20 accessions / species (Hélène, CIAT data)
    - PCA: separation of the 4 species, overlap of the phenotype between African/Asian and separation wild/cultivated
    - phenotyping IR64/Nipponbare/W1654/Tog5681/B88 (Hélène's data from IRD)
* **Supplementary figure/table 1**:
    - Timing of developmental transitions (did not detect heterochrony)
    - How does the lack of heterochrony relate to gene expression?

### 2. General RNAseq Results

* **Supplementary figure/table 2**:
    - General RNAseq stuff:
        + Mapping QC, to justify mapping against *O. sativa japonica* (mapping stats for all libraries against OS, mapping stats for each library against draft genome)
    - "Validation" of stages by qPCR
* **Figure 2**:
    - Sample clustering & heatmap
    - species cluster together, stages cluster within species
    - african species cluster together, asian species cluster together -> domestication changes have been small compared to species-level changes
* **Supplementary figure/table 3**:
    - Differential expression between species (either by comparing stages between species or comparing species)
    - these results are not necessarily related to branching
    - There are 1000s of DE genes between species, but this analysis is confounded by mapping bias against the O. sativa japonica reference
    - Can we detect coding differences between genes from the RNAseq data? (Tom has done this analysis)
* **Supplementary figure/tables 4**:
    - DE analysis, all genes
    - use one of these comparisons (e.g. stage:domestication interaction) to justify focus on TF families

### 3. Evolution of transcription factor expression

* We have to decide how to discuss this. Do we present individual genes, or talk on the family level?
* **Figure 3**:
    - TFs differentially expressed between stages in all species.
    - These are the "basal" TFs involved in branch meristem maturation in Oryza spp.
* **Figure 4**:
    - Interaction between continent and stage (TFs where the DE between stages is different for African and Asian rice).
    - Discuss these TFs wrt. phenotypic differences between African and Asian rice that could be related to differential expression of these genes. African species have fewer secondary branches.
* **Figure 5**:
    - Each domesticated accession wrt. to its wild relative (interaction between stage and accession). The wild accessions have fewer primary branches compared to the domesticated accessions.
        + rufipogon vs. indica
        + barthii vs. glaberrima
        + rufipogon vs. japonica (*maybe*)
        + rufipogon vs. indica + japonica (*maybe*)
    - These TFs are putative targets of selection in domestication of indica and glaberrima (+/- japonica)
* **Figure 6**:
    - Parallel evolution: interaction between stage and domestication (+/- japonica)
    - These are essentially the TFs from **figure 3** where the DE between stages is dependent on whether the species is domesticated or wild
    - This analysis treats rufipogon and barthii as a wild "pool" and indica and glaberrima (+/- japonica) as the domesticated "pool". Logically, this does not make sense because rufipogon and barthii are not the same. However, it's statistically valid, and it gives us some interesting genes.
    - We only have two domestications to compare. Anything we pick out of this analysis will be "associated" with domestication, but we can't make functional / mechanistic claims. This goes for all the comparisons in our paper, e.g. rufipogon vs. glaberrima we are only comparing one accession of each, we don't know what's happening in all the other glaberrima accessions.

### 4. qPCR screen of TFs

* Pick 1/2 TF families from **3** and test their expression over the whole range of samples (RM, BM, SM, FM). AP2? MADS?
* **Figure 7**:
    - qPCR results from family X
    - which families do we have qPCR results for?
* **Supplementary figure/table 5**:
    - qPCR of certain candidate genes identified in LMD paper
    - do we have this data already?

### 5. Something general

* Would be good to finish the results/discussion with something general. These are ideas, please add more:
    - Can we see any TF family-level effects? z-score / enrichment, etc.
    - Is there anything interesting about redundancy? Different subsets of genes from the same family being selected in the different domestications? Tree of genes from family x and the selection response in different accessions?

##  Conclusion

Domestication / branching complexity correlates with AP2 gene expression

## To discuss (not necessarily in the paper)

* Which journal

##

<div id="refs"></div>
