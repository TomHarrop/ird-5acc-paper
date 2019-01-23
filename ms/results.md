## Results

### Parallel changes in panicle architecture between wild and domesticated accessions

To measure the diversity of panicle architecture, we phenotyped 91 rice accessions (Table S1 **Plantinfo**), including wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*), using `P-TRAP` image analysis software for automated measurement of traits [Figure 1A **phenotype-pca**; Figure S1 **PanicleStructure**; Table S1 **Plantinfo**; Table S3 **PanicleTraitsPhenotyping**; @al-tamPTRAPPanicleTrait2013].
Principal components analysis (PCA) of the phenotyping data identified a major coordinate (PC1) that accounts for 46.5% of variability (Figure 1B **phenotype-pca**).
PC1 separates domesticated and wild accessions, but not Asian and African accessions, and is the only component that separates panicles from different accessions (Figure S2 **phenotype-pca-all-pc**).
Spikelet number, secondary branch number and primary branch number have the highest loadings on PC1, whereas length traits have lower absolute loading on PC1 (Figure 1C **phenotype-pca**).
For all species, spikelet number correlates more with secondary branch number than it does with primary branch number.
Primary branch number correlates with spikelet number more in wild than domesticated species, but this correlation is weaker in Asian species than in African species.
Primary and secondary branch numbers do not correlate, suggesting they are controlled by different genetic mechanisms (Figure 1D **phenotype-pca**).
Our phenotypic analysis indicates similar changes in panicle architecture between wild and domesticated accessions in the independent African and Asian domestication processes, which presumably result from parallel, artificial selection on panicle architecture.
Spikelet number, secondary branch number and primary branch number are the main contributors to these differences in panicle architecture, and these phenotypes are all related to axillary meristem formation and fate transition 
[@teoNewInsightsRegulation2014; @zhangMolecularControlGrass2014].

### Measurement of gene expression in developing panicles

We investigated gene expression changes associated with the diversity of panicle architecture and differences between the Asian and African domestication processes via RNA sequencing (RNAseq).
We used a single accession each of domesticated Asian rice (*O sativa indica* IR64) and its wild relative (*O. rufipogon* W1654), and domesticated African rice (*O. glaberrima* Tog5681) and its wild relative (*O. barthii* B88).
We also included *O. sativa japonica* cv. Nipponbare as the genomic reference accession.
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S3 **phenotype-all-varieties**).
To measure whole-transcriptome gene expression in these accessions, we collected immature panicles from each accession at four developmental stages:
rachis meristem (RM);
indeterminate meristem (IM), including panicles displaying primary branch initiation, elongation of primary branches and axillary meristem initiation;
determinate meristem (DM), including panicles wherein axillary meristems had transitioned into early spikelet differentiation;
and floret meristem (FM), with early differentiation of floral organs
(Figure S4a **qpcr-confirms-sampling**).
We first confirmed staging of the panicles by extracting RNA from pooled immature panicles at each stage and measuring expression of markers of panicle development by quantitative real-time RT-PCR (qPCR) (Figure S4b **qpcr-confirms-sampling**).
Because branching complexity is related to branch meristem establishment and meristem fate transition [@kyozukaControlGrassInflorescence2014], and secondary branch number and spikelet number contribute to differences between wild and domesticated accessions (Figure 1), we chose the IM and DM stages for RNAseq.
cDNA libraries for sequencing were constructed from rRNA-depleted RNA samples from three biological replicates at both stages for all five accessions.
After mapping the reads against the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007], we obtained an average of more than 20 million uniquely mapped reads within exons for each accession, including African rice species (Table S4 **mapping-statistics**).
Pairwise distances between samples, calculated on the number of reads per gene from all detected genes, grouped samples first by stage, then by accession, and then by continent (Figure S5 **distance-heatmap**).
We did not observe grouping by domestication status, suggesting that transcriptome-wide changes during domestication are subtle compared to differences between species.

### *AP2/EREBP*-like transcription factors are core regulators of panicle branching

We used PCA on transformed raw counts to investigate general patterns of variation in the transcriptomes that may be associated with differences in panicle architecture. The first four PCs split different combinations of rice species (Figure 2 **transcriptome-pca**).
PC1–PC4 may relate to species-specific differences unrelated to panicle architecture, or mapping biases introduced by mapping all libraries against the *O. sativa japonica* reference.
In contrast, PC5 separates developmental stages across all five accessions, although separation is weaker in *O. sativa indica* (Figure 2 **transcriptome-pca**). 
*APETALA2* and ethylene-responsive element binding protein (*AP2/EREBP*)-like genes and MADS-box genes are enriched at the extremes of PC5 (*p*~adj~ = 0.004 for both, GSEA permutation *t*-test; Table S5 **PC5-TF-enrichment-gsea**).
Generally, *AP2/EREBP*-like genes that contribute to PC5 are more highly expressed at the IM stage, and MADS-box genes are more highly expressed at the DM stage (Figure 3 **HB-AP2-heatmap**).
Most of the *AP2/EREBP*-like genes that are highly expressed at the IM stage belong to the ERF clade.
*NAM*, *ATAF1* and *CUC2*-domain (NAC) transcription factor genes are also more highly expressed at the IM stage, and *SQUAMOSA* PROMOTER BINDING PROTEIN-like (SBP) genes including *WFP* and *SPL7* are more highly expressed at the DM stage [@miuraOsSPL14PromotesPanicle2010; @wangCoordinatedRegulationVegetative2015].
There is a group of homeobox (HB) genes that is highly expressed in the DM samples, and a larger group that is highly expressed at the IM stage, and ten out of twenty of these genes belong to the HD-ZIP clade [Figure S6 **NAC-HB-SPL-heatmap**; @jainGenomewideIdentificationClassification2008].
Co-regulation of members of transcription factor (TF) families, sometimes at the clade level, highlights the redundant or overlapping functions of TF families in meristem establishment and fate transition.

To identify the core set of genes involved in axillary meristem determination in all five accessions, we used differential expression (DE) tests to find genes that were up- or down-regulated between stages across all accessions.
Positive log~2~-fold change values (L~2~FCs) indicate higher expression in DM than in IM.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1, including 115 genes up-regulated in DM samples and 78 genes down-regulated in DM samples (Table S6 **DE-genes-stages**).
There was an enrichment of transcription factor (TF) genes in the list of 193 DE genes (38 TF genes; *p* = 1.9 × 10^-9^, hypergeometric test), including *FZP*, *LHS1*, *LAX1*, *PAP2*, and *MFO1*, which regulate inflorescence architecture in rice [@komatsuLAX1FRIZZYPANICLE2001; @baiRegulatoryRoleFZP2016; @baiDuplicationUpstreamSilencer2017; @kobayashiPANICLEPHYTOMER2PAP22010; @kobayashiInflorescenceMeristemIdentity2012; @khandayRiceLHS1OsMADS12013; @ohmoriMOSAICFLORALORGANS12009a; @huangVariationRegulatoryRegion2018].
This indicates that RNAseq of developing panicles at the IM and DM stage identifies genes that control branching.
The list of 193 DE genes included 9 MADS-box genes and 10 *AP2/EREBP*-like genes (*p*~adj~ 5.2 × 10^-11^ and 1.4 × 10^-8^ respectively, hypergeometric test; Table S6 **DE-genes-stages**), in agreement with the prominence of those families in PC5.

The DE and PCA results are consistent with the role of transcriptional regulation in panicle branching, and highlight a set of candidate core regulators of axillary meristem determination and branching that are conserved in rice.
The pattern of expression of *AP2/EREBP*-like genes may indicate a role in the promotion of indeterminate axillary meristem identity or suppression of the transition from axillary meristem to spikelet meristem.
MADS-box and HB genes may have an inverse role as promoters of determinate meristem.

To test the role of *AP2/EREBP*-like genes in the control of panicle architecture, we phenotyped panicles from two loss-of-function mutants.
The *crl5* mutant of the *AP2*-like gene *PLETHORA 8* [*PLT8*; @kitomiAuxinResponsiveAP22011] produced panicles with a shorter rachis with fewer primary branches (Figure 4 **panicle-mutants**; Table S7 **PanicletraitsPhenotypingAP2**), consistent with a peak of *PLT8* expression in rachis meristem tissues from *O. sativa japonica* [Figure S7 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
Panicles from the *smos1* mutant of *ERF142* [@ayaNovelAP2TypeTranscription2014] have a reduced number of primary and secondary branches, and fewer spikelets (Figure 4 **panicle-mutants**; Table S7 **PanicletraitsPhenotypingAP2**).
*ERF142* expression is highest in primary branch and elongating primary branch meristem tissues in *O. sativa japonica* [Figure S7 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
Although neither gene was differentially expressed in our RNAseq dataset, these phenotypes support the involvement of *AP2/EREBP*-like genes in control of panicle architecture.

### *AP2/EREBP*-like gene expression is associated with panicle diversity and domestication

To identify common patterns of expression of transcriptional regulators related to variation in panicle phenotypes, we used soft clustering of scaled L~2~FCs between IM and DM using the subset of annotated TF genes that were detected in our RNAseq dataset.
We did not include the reference accession, *O. sativa japonica* cv. Nipponbare, because it is not a high-yielding cultivar compared to other Asian rice accessions (Figure S3 **phenotype-all-varieties**).
We recovered six clusters comprised of a total of 119 genes (Figure 5 **cluster-phenotype-corr**; Table S8 **clustered-genes**).
We calculated correlations between the mean L~2~FC value of genes in each cluster and PC1 in the phenotyping data, and between mean L~2~FC and the number of secondary branches and spikelets from repeat panicle phenotyping for the accessions used for RNAseq (Figure 5b **cluster-phenotype-corr**; Figure 1 **phenotype-pca**; Figure S8 **phenotyping-mpl**).

Clusters 3, 4 and 6 correlate with spikelet number and secondary branch number, but not with PC1, meaning that the L~2~FC of genes in those clusters does not correlate with the phenotypic differences between wild and domesticated accessions.
Clusters of genes with high L~2~FC in *O. sativa indica* have a positive correlation with SBN and SpN, whereas clusters of genes with low L~2~FC in *O. sativa indica* have a negative correlation with SBN and SpN.
Cluster 4 had an enrichment of HB genes (9 out of 31 genes; *p*~adj~ = 2.5 × 10^-4^).
It also contained three MIKC^C^-type MADS-box genes (*LHS1*, *MFO1* and *MADS14*), which promote spikelet meristem determination [@jeonLeafyHullSterile12000; @ohmoriMOSAICFLORALORGANS12009a; @kobayashiInflorescenceMeristemIdentity2012], and three *AP2/EREBP*-like genes including *INDETERMINATE SPIKELET 1* (*IDS1*), which also controls inflorescence architecture [@chuckFloralMeristemInitiation2008; @leeTwoAP2Family2012]. 
L~2~FCs of genes in this cluster are low in *O. sativa indica*, high in *O. glaberrima* and intermediate in the two wild accessions.
Although these genes may be involved in regulation of panicle complexity, their expression did not appear to have changed in parallel in the two domestications.
L~2~FCs of genes in clusters 3 and 6 change between accessions from the two continents.
In cluster 3, L~2~FCs are higher in African species than in Asian species, meaning that the genes are more highly expressed in DM stages in African species.
Genes in cluster 6 have the opposite pattern, with higher L~2~FCs in Asian species compared to African species.
Cluster 3 contains the *LAX PANICLE 1* (*LAX1*) and *FLO-LFY HOMOLOG OF RICE* (*RFL*), which are involved in axillary meristem establishment and outgrowth and promotion of indeterminate meristematic activity in rice respectively [@komatsuLAX1FRIZZYPANICLE2001; @ikeda-kawakatsuABERRANTPANICLEORGANIZATION2012].
Their higher expression at the DM stage in panicles from both wild and domsticated African accessions could be associated with a reduced number of spikelets.

Clusters 1, 2 and 5 correlated with the main principal component (PC1) in the phenotyping data, which separates wild and domesticated species independent of continent.
Clusters 1 and 5 are also positively correlated with spikelet number and secondary branch number, whereas cluster 2 is negatively correlated. 
The correlation of L~2~FC patterns with PC1 suggests that genes in these clusters may be associated with changes in panicle architecture between wild and domesticated species.
Cluster 1 and cluster 5 both had a positive correlation with PC1.
L~2~FCs are higher in domesticated accessions for genes in cluster 1, meaning that they are more highly expressed at the DM stage in domesticated accessions.
Genes in cluster 2 have lower L~2~FCs in domesticated accessions, meaning that they are more highly expressed at the DM stage in wild accessions.
Cluster 2 also had the strongest negative correlation with PC1, and low L~2~FCs in *O. sativa indica*. 
The lower L~2~FCs in domesticated accessions could implicate cluster 2 genes in promotion of determinate meristem fate, because their lower expression at the DM stage may result in more activity of indeterminate axillary meristems.
Eight *GROWTH-REGULATING FACTOR1* (GRF) family genes, which are involved in the regulation of cell proliferation [@kimRegulationPlantGrowth2015], were detected in our dataset, and three of them were present in cluster 2 (*p*~adj~ = 7.6 × 10^-3^, hypergeometric test).
In contrast to cluster 2, most of the genes in cluster 5 have L~2~FCs close to zero in *O. sativa indica*, and negative L~2~FCs in the other accessions (Figure S9 **cluster-5-details**), meaning that the expression of these genes decreases at the DM stage in all accessions except *O. sativa indica*.
The lack of repression of cluster 5 genes and to a lesser extent cluster 1 genes at the DM stage in *O. sativa indica* could result in more branching via the promotion of indeterminate/axillary meristem identity.
Cluster 5, which had the highest mean L~2~FC in *O. sativa indica*, had an enrichment of *AP2/EREBP*-like genes (6 out of 17 genes; *p*~adj~ = 7.6 × 10^-3^, hypergeometric test), and cluster 1 also contains 3 *AP2/EREBP*-like genes (Table S8 **clustered-genes**).
We used qPCR to confirm these patterns in all four stages of each accession for all *AP2/EREBP*-like genes in cluster 5 (Figure S10 **fluidigm-ap2-hb**).
The enrichment of *AP2/EREBP*-like genes in cluster 5 and the presence of 3 *AP2/EREBP*-like genes in cluster 1 suggests that their pattern of expression is associated with differences in panicle architecture across wild and domesticated accessions.

To find TF genes associated with parallel changes in panicle architecture during domestication, we tested the stage × domestication interaction for *O. rufipogon*, *O. sativa indica*, *O. barthii* and *O. glaberrima* at an FDR of 0.1 (Table S9 **DE-genes-interaction**).
We detected 19 genes with a stage × domestication interaction, including nine *AP2/EREBP*-like genes (*p* = 4.4 × 10^-7^, hypergeometric test; Figure 6A **dom-genes-plot**).
Seven of these *AP2/EREBP*-like genes are in the top 10% of all genes by absolute loading on PC5, with higher expression in IM than DM (Figure 3 **MADS-AP2-heatmap**), and higher L~2~FC in wild than domesticated species (Figure 6A **dom-genes-plot**).
These genes are putative targets of parallel selection on panicle architecture that occurred during domestication.
*AP2/EREBP*-like genes were also prominent when we tested the stage × accession interaction separately for each domestication (12 out of 85 genes in Asian accessions; 8 out of 50 genes in African accessions; Table S9 **DE-genes-interaction**).
Consistent with the presence of *IDS1* in cluster 4, it was also was DE in both Asian and African domestications, although the direction of change was different (Figure 6B **dom-genes-plot**).
Genes with this pattern of expression in the four accessions may have also been targets of selection on panicle architecture, but evolved divergently.

The prominence of *AP2/EREBP*-like genes among putative core regulators of branching in all four *Oryza* species, and among genes associated with differences between wild and cultivated accessions, suggest that they were key targets of artificial selection for improvement in panicle architecture, and were involved in changes to the regulatory network controlling branching that occurred during domestication.
