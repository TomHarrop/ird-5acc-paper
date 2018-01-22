
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

### 1. Phenotypic analyses of African and Asian domesticated and wild rice

* Parallel evolution of phenotype during African and Asian domestications
* Correlations between *e.g.* spikelet number and branching complexity
    - these correlations justify the choice of PBM and SM for sequencing
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
* **Figure 2**:
    - Sample clustering & heatmap
    - species cluster together, stages cluster within species
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
    - Discuss these TFs wrt. phenotypic differences between African and Asian rice that could be related to differential expression of these genes
    - I'm not sure about this section, maybe supplementary data, or remove
* **Figure 5**:
    - Each domesticated accession wrt. to its wild relative (interaction between stage and accession)
        + rufipogon vs. indica
        + barthii vs. glaberrima
        + rufipogon vs. japonica (*maybe*)
        + rufipogon vs. indica + japonica (*maybe*)
    - These TFs are putative targets of selection in domestication of indica and glaberrima (+/- japonica)
* **Figure 6**:
    - Parallel evolution: interaction between stage and domestication (+/- japonica)
    - These are essentially the TFs from **figure 3** where the DE between stages is dependent on whether the species is domesticated or wild
    - This analysis treats rufipogon and barthii as a wild "pool" and indica and glaberrima (+/- japonica) as the domesticated "pool". Logically, this does not make sense because rufipogon and barthii are not the same. However, it's statistically valid, and it gives us some interesting genes.

### 4. qPCR screen of TFs

* Pick 1/2 families from **3** and test their expression over the whole range of samples (RM, BM, SM, FM). AP2? MADS?
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

Reshaping of the ~maize~ rice transcriptome by domestication :P

## To discuss (not necessarily in the paper)

* Which journal