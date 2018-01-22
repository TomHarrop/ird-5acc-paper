
## Introduction

* evo devo and gene expression
* domestication and gene expression (maize, tomato transcriptomes)
* details about rice domestication
* **aims**:
    - changes in phenotype, gene expression and metabolism caused by / associated with domestication
    - common genes that respond to selection in different domestications
    - genes underlying the branching process in oryza spp.

## Results/discussion

### 1. Phenotypic analyses of African and Asian domesticated and wild rice

* Parallel evolution of phenotype during African and Asian domestications
* Correlations between *e.g.* spikelet number and branching complexity
    - because of these correlations we can justify PBM and SM
* Figure 1:
    - phenotyping 20 accessions / species (Hélène, CIAT data)
    - PCA: separation of the 4 species, overlap of the phenotype between African/Asian and separation wild/cultivated
    - phenotyping IR64/Nipponbare/W1654/Tog5681/B88 (Hélène's data from IRD)
* Timing of developmental transitions:
    - did not detect heterochrony (Supp Data/ Supp Figure)
    - therefore investigate gene expression

### 2. Gene expression 1: qPCR screen

* Figure 2 / table 1: qPCR of certain candidate genes identified in LMD paper
    - Which genes? Highlight some genes DE between stages PBM and SM to justify focus on these stages for RNAseq

### 3. Gene expression 2: RNAseq Results

* General RNAseq stuff (in supp data):
    - qc results figure: Justify mapping against *O. sativa japonica* (read mapping stats, for each genome?)
    - Sample clustering/heatmap: species cluster together, stages cluster within species (compare with metabolite clustering for African species)
* Figure 5: RNAseq analysis of metabolite pathways. Heatmap / something else showing expression differences of genes involved in one/two of the pathways we discussed in 3.
* Figure 6 (or supp data): Differential expression between species (either by comparing stages between species or just comparing species)
    - we're ignoring branching, but maybe that fits with the metabolic data
    - Link this to metabolites by focussing on differences in genes from a metabolic pathway.
    - Can we detect coding differences between genes in these pathways from the RNAseq data? (n.b. we've done this experiment, I just have to find the results...)
    - Too many genes to go into a lot of detail in this section.
    - RNAseq mapping also problematic here
* Figure 7: Genes differentially expressed between stages in all species.
    - These are the "basal" genes involved in branch meristem maturation in Oryza.
* Figure 8: Interaction between continent and stage.
    - Genes where the DE between stages is different for African and Asian rice.
    - What are the phenotypic differences between African and Asian rice that could be related to differential expression of these genes?
    - Not sure about this section, maybe supplementary data, or remove completely?
* Figure 9: Each domestication wrt. to its wild relative (interaction between stage and species)
    - These are genes that may have been selected in domestication of japonica, indica and glaberrima
    - Any or all of these comparisons:
        + rufipogon vs. japonica
        + rufipogon vs. indica
        + rufipogon vs. (indica + japonica)
        + barthii vs. glaberrima
* Figure 10: Parallel evolution: interaction between stage and domestication
    - These are the genes where the DE between stages is dependent on the species is domesticated or wild, treating rufipogon and barthii as a wild "pool" and japonica, indica and glaberrima as the derived "pool"
    - n.b. logically this is weird, because rufipogon and barthii are not the same (Yves's point), but statistically it is OK, and it gives us some interesting genes.

### AP2 genes?

##  Conclusion

?

## Points à discuter (not necessarily in the paper)

* Which journal

### Metabolites

* Metabolite results suggest that the Asian species are more divergent than the African species. Is this because the domestication of African rice happened more recently, so there has been less phenotypic evolution? How would the effects of evolution be different between gene regulation and metabolic pathways? Evo devo suggests that regulatory changes produce morphological change, but how does this apply to metabolism?

* In the African species' transcriptomes, the clustering is **within** species, so this might support the 'faster evolution of transcriptomes than metabolomes' idea. It might also suggest that the metabolic differences do not relate to differences in inflorescence development *per se*, because they seem to be more about organismal differences than stage-specific differences. 

### RNAseq

* How can we present "sets" of genes? Some Wald tests result in hundreds of DE genes.
    - option 1: filter aggressively (*p*~adj~, L~2~FC threshold etc.)
    - option 2: some sort of geneset-level overview like KEGG pathways, pfam...?
    - other options?
