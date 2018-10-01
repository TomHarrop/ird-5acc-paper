## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice, followed by whole-transcriptome sequencing (RNASeq) [and analysis of single-gene mutants with defects in panicle branching]. **expand slightly**

### Panicles from domesticated accessions produce more primary branches, secondary branches and spikelets

To measure the diversity of panicle architecture, we phenotyped **n** accessions of wild Asian rice (*Oryza rufipogon*; Figure 1A), domesticated Asian rice (*Oryza sativa*; Figure 1B), wild African rice (*Oryza barthii*; figure 1C) and domesticated African rice (*Oryza glaberrima*; Figure 1D), using `P-TRAP` for automated measurement of traits [@al-tamPTRAPPanicleTrait2013]. The first principal component in the phenotyping data accounts for 46.5% of variability, and separates domesticated and wild accessions but not Asian and African accessions (Figure 1E). Spikelet number, secondary branch number and primary branch number have the highest loadings on PC1, suggesting that these factors are the main differences between panicles from wild and domesticated accessions. Lower-order components do not separate panicles from different accessions (Figure S1).

### Differences in gene expression profiles between wild and domesticated panicles

To investigate gene expression associated with variation in panicle architecture, we chose a single accession for each of *O. rufipogon*, *O. sativa* japonica (Nipponbare), *O. sativa* indica (IR84), *O. barthii* (**?**) and *O. glaberrima* (**TOG...**) for RNA sequencing (RNAseq).

To test whether these accessions are consistent with species-wide patterns of panicle architecture, we performed further panicle phenotyping on **n** panicles from **m** plants from each accession. The three domesticated accessions produce more spikelets than their wild relatives (Figure S2), and there was a stronger correlation between secondary branch number and spikelet number in the domesticated accessions than in the wild accessions, suggesting **...** (Figure S3). **Model of panicles goes here.**

To measure whole-transcriptome gene expression, we collected developing panicles from each accession at rachis meristem (RM), branch meristem (BM), spikelet meristem (SM) and **floret?** meristem (FM) stages (Figure S4a). We extracted RNA from all stages, and measured expression of *LHS1* and *G1L5* (*TAWAWA*) with quantitative RT-PCR (qRT-PCR) to confirm staging (Figure S4b). We used RNA from the BM and SM stages for rRNA depletion and construction of cDNA libraries for RNAseq. We obtained an average of **x** million 100b read pairs per sample (Table S1), which we mapped against the MSU Release 7.0 of the annotation of the *O. sativa* japonica cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007]. Unique read mapping percentages were above **X**% even for African rice species (Table S1). Pairwise distances between libraries were lower for samples from the same stage, accession and continent, in that order (Figure S5), indicating that transcriptomic changes during domestication are minor compared to transcriptome differences between species. We used differential expression (DE) analysis to identify genes that were up- or down-regulated between stages across all accessions. 193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate cutoff 0.1 (Table S2). This list included **genes x, y and z** [e.g. pick from *FZP*, *LHS1*, *LAX1*, *PAP2*, *MFO1*, ...], which control inflorescence architecture in rice [**add citations for chosen genes**]. 

## APETALA2 and ethylene-responsive element binding proteins transcription factors are differentially expressed between stages

We used PCA to investigate the patterns of variation in the transcriptomes. The first four PCs split different combinations of rice species (Figure 2). These components may relate to species-specific differences unrelated to panicle architecture, or mapping biases caused by mapping all libraries against the *O. sativa* japonica reference. In contrast, PC5 separates developmental stages across all five accessions except *O. sativa* indica (Figure 2). The family of APETALA2 and ethylene-responsive element binding protein (AP2/EREBP)-like genes are enriched at the extremes of PC5, suggesting that they contribute to differences between BM and SM (Figure 3). Homeobox, MADS, NAC and SBP genes are also enriched in PC5 (Figure 3; Figure S5). 

## Figure legends

**Figure 1**. Panicles from domesticated accessions produce more branches and spikelets than panicles from wild accessions. We used spread panicles from **A** *O. rufipogon*, **B** *O. sativa*, **C** *O. barthii* and **D** *O. glaberrima* to measure panicle phenotypes with P-TRAP [@al-tamPTRAPPanicleTrait2013]. The first principal component (PC1) in the panicle phenotye data accounts for 46.5% of variability and separates wild and domesticated accessions (**E**), and spikelet number (SpN), secondary branch number (SBN) and primary branch number (PBN) had the highest loadings on PC1 ( **F**). RL: Rachis length; PBL: Primary branch length; PBIL: Primary branch internode length; SBL: Secondary branch length; SBIL: Secondary branch internode length; TBN: Tertiary branch number; PL: Panicle length.

**Figure 2**. Principal component 5 (PC5) separates libraries by developmental stage, and explains 5.4% of total variability. The first four components explain 51.7% of variability, and separate libraries by species.

**Figure 3**. AP2/EREBP and homeobox (HB) transcription factors change expression between BM and SM. **A**. AP2/EREBP and HB genes are distributed at the extremes of genes ranked on PC5 (GSEA test *p*~adj~ 0.004 and 0.021, respectively **which test? hypergeometric?**). We retained the highest 10% of genes by absolute value on PC5 (**check wording**). **B**. Most AP2/EREBP genes that pass the cutoff are more highly expressed (**what's on the scale bar?**) in the BM. Three of the four AP2/EREBP genes that are more highly expressed in the SM belong to the AP2 subfamily. Genes that are more highly expressed in BM mainly belong to RAV, DREB and ERF subfamilies. **C**. Most HB genes that pass the cutoff are more highly expressed in the SM. Ten out of twenty of these genes belong to the HD−ZIP IV subfamily (**are there stats for this?**).

**Figure S1**. Principal component analysis (PCA) of panicle phenotyping data showing components 1–4. PC1 accounts for 46.5% of variability and separates panicles from domesticated and wild accessions. The lower ordinates do not separate panicles by species.

**Figure S2**. Detailed phenotyping of five Oryza accessions. The three domesticated accessions produce more spikelets than their wild relatives. In comparison to *O. sativa* japonica, *O sativa* indica produces more secondary branches.

**Figure S3** Correlation between the main traits that determine panicle phenotype. A. Primary Branch Number and Spikelet Number correlates only slightly mainly in wild species. B. Secondary Branch Number and Spikelet number
highly correlates in cultivated species, less in wild species. C.
Primary and Secondary Branches poorly correlates, suggesting that they
are controlled by different genetic mechanisms

**Figure S4**. **A** Photos of panicle samples. **B** *LHS1* / *G1L5* qRT-PCR.

**Figure S5**. MADS and SBP genes are more highly expressed in the SM. NAC genes are more highly expressed in the BM.

## Tables

**Table S1**. Read and mapping statistics for all RNAseq libraries.

**Table S2**. Genes diffentially expressed between stages across all species.