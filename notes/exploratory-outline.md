
## Introduction
blababababababab
* evo devo and gene expression
* transcription factors / evolution of gene regulation
* domestication and gene expression (maize, tomato transcriptomes)
* details about rice domestication
* branching development in rice
* **aims**:
    - changes in phenotype, gene expression caused by / associated with domestication
    - common genes that respond to selection in different domestications
    - evolution of transcription factor family expression
    - does selection for the same phenotype affect expression of the same genes / sets of genes

## Results/discussion

**EVIDENCE**

### 1. Phenotypic analyses of African and Asian domesticated and wild rice

* Parallel evolution of phenotype during African and Asian domestications
* breeders selected for similar panicle phenotype in the 2 independent domestications
* Correlations between *e.g.* spikelet number and branching complexity
    - these correlations justify the choice of PBM and SM for sequencing
    - also, "evolutionary diversity in Solanaceae inflorescence complexity is determined by subtle modifications of transcriptional programs during a critical transitional window of meristem maturation" [@Lemmon_Evolutioninflorescencediversity_2016]
* **Figure 1**:
    - photo of mature panicles
    - phenotyping 20 accessions / species (Hélène, CIAT data)
    - PCA: separation of the 4 species, overlap of the phenotype between African/Asian and separation wild/cultivated
    - phenotyping IR64/Nipponbare/W1654/Tog5681/B88 (Hélène's data from IRD)
* **Supplementary figure/table 1**:
    - Timing of developmental transitions (did not detect heterochrony)
    - How does the lack of heterochrony relate to gene expression?

### 2. General RNAseq Results

* **Supplementary figure/table 2 (exclude this)**:
    - General RNAseq stuff:
        + cross species mapping is acceptable
        + Mapping QC, to justify mapping against *O. sativa japonica* (mapping stats for all libraries against OS, mapping stats for each library against draft genome)
    - "Validation" of stages by qPCR
* **Supplementary figure/tables 3**:
    - DE analysis, all genes
    - use one of these comparisons (e.g. stage:domestication interaction) to justify focus on TF families
* **Figure 2**:
    - Sample clustering & heatmap
    - species cluster together, stages cluster within species
    - african species cluster together, asian species cluster together -> domestication changes have been small compared to species-level changes
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
        + rufipogon vs. japonica (*exclude*)
        + rufipogon vs. indica + japonica (*exclude*)
    - These TFs are putative targets of selection in domestication of indica and glaberrima
    - genes with a +ve L2FC for the stage:species interaction fit our hypothesis (if the wild species is lower level in the species factor)
    - there is some non-overlap of glab / indica genes. Of course selection could have acted on different genes in different domestications
* **Figure 6**:
    - Parallel evolution: interaction between stage and domestication (+/- japonica)
    - These are essentially the TFs from **figure 3** where the DE between stages is dependent on whether the species is domesticated or wild
    - This analysis treats rufipogon and barthii as a wild "pool" and indica and glaberrima (+/- japonica) as the domesticated "pool". Logically, this does not make sense because rufipogon and barthii are not the same. However, it's statistically valid, and it gives us some interesting genes.
    - We only have two domestications to compare. Anything we pick out of this analysis will be "associated" with domestication, but we can't make functional / mechanistic claims. This goes for all the comparisons in our paper, e.g. rufipogon vs. glaberrima we are only comparing one accession of each, we don't know what's happening in all the other glaberrima accessions.
* **Supplementary figure/table 4**:
    - Differential expression between species (either by comparing stages between species or comparing species)
    - these results are not necessarily related to branching
    - There are 1000s of DE genes between species, but this analysis is confounded by mapping bias against the O. sativa japonica reference
    - Can we detect coding differences between genes from the RNAseq data? (Tom has done this analysis)
* **Supplementary figure/table 5**:
    - qPCR results from interesting TFs
    - done for AP2, ALOG, SPL, MADS (SEP, SVP, AP1 subfamily)
* **Supplementary figure/table 6**:
    - qPCR of certain candidate genes identified in LMD paper
    - done for SPL, ALOG and some genes related to hormones

### 3. Correlation between phenotype and gene expression

* **Figure 7**
    - correlation between PCs from the phenotypic data and L2FCs from the RNAseq data
    - genes where the correlation between PCx and L2FC is +ve would suggest higher expression in the SM results in a higher score on PCx
    - not sure which genes or which component yet

### 4. Patterns of gene expression in the RNAseq data

* **Figure 8**
    - PCA of RNAseq alone (no phenotype),
    - Most of the differences are explained by the inter-species
    - components 1 to 4th explain interspecies diff, component 5 explains stage differences, the rest is noise
    - PC5:
        + This component splits the two stages and it is mostly driven by the wild species, while nipponbare and especially indica are much closer to 0.
        + The wild species display a much stronger burst of MADS-box related to differientiation
        + The AP2 story might be more complicated, but indeed many of them are strongly down-regulated in most species, but they behave differently in indica, nipponbare

* **Figure 9**
    - clustering of TF genes
    - regulatory "modules" that have co-evolved within and between species, e.g. cluster 6 contains genes that have a higher L2FC in indica and glaberrima than in barthii and rufipogon
    - no family survives bonferroni correction in cluster 6
    - AP2 genes are enriched in cluster with high +ve L2FC in indica.
    - Hélène is currently looking for biologically interesting stuff in the clusters

### 5. Evolution of CRMs

- align / analyse a few genes
- link the DE genes with SNP variations in the different species in the promoter
- Based on IRRIGIN dataset (O. barthi, O. glaberrima, O. rufipogon and O. sativa)
- **choose genes**

**CONCLUSION/HYPOTHESIS**

- there is a group of genes that are repressed at the SM stage in wild species (-ve L2FC PBM vs. SM). These genes "promote indeterminate state". These genes are *not* repressed in the domesticated species, therefore the indeterminate state is maintained for longer
    + **domestication/branching complexity is [partially] controlled by repression of genes that promote indeterminate state at the SM stage**
- it may not be the same genes in glaberrima and indica, but it is the general mechanism that is interesting
- genes that show the same pattern in both domesticated species -> parallel evolution
- family-level comparisons are interesting too but not the whole story, e.g. an AP2 TF could regulate a MADS TF etc.
-  is there anything interesting about redundancy? Different subsets of genes from the same family being selected in the different domestications? Tree of genes from family x and the selection response in different accessions?

##

<div id="refs"></div>
