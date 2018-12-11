## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice. We selected five accessions for whole-transcriptome sequencing (RNAseq). Our analysis reveals an association between expression of *AP2/EREBP*-like genes and panicle architecture.

### Parallel changes in panicle architecture between wild and domesticated accessions

To measure the diversity of panicle architecture, we phenotyped 91 accessions of wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*), using `P-TRAP` image analysis software for automated measurement of traits [Figure 1A **phenotype-pca**; Figure S1 **PanicleStructure**; Table S3 **PanicleTraitsPhenotyping**; @al-tamPTRAPPanicleTrait2013].
The first principal component (PC1) in the phenotyping data accounts for 46.5% of variability.
PC1 separates domesticated and wild accessions, but not Asian and African accessions, and spikelet number, and is the only component that separates panicles from different accessions (Figure S2 **phenotype-pca-all-pc**).
Secondary branch number and primary branch number have the highest loadings on PC1 (Figure 1B–C **phenotype-pca**).
For all species, spikelet number correlates more with secondary branch number than it does with primary branch number.
Primary branch number correlates with spikelet number more in wild than domesticated species, but this correlation is weaker in Asian species than in African species.
Primary and secondary branch number do not correlate strongly, suggesting they are controlled by different genetic mechanisms (Figure 1D **phenotype-pca**).
Our phenotypic analysis indicates similar changes in panicle architecture between wild and domesticated accessions in the independent African and Asian domestication processes, which presumably result from parallel, artificial selection on panicle architecture.
Spikelet number, secondary branch number and primary branch number are the main contributors to differences in panicle architecture.
These phenotypes are all related to axillary meristem formation and fate transition.

### Measurement of gene expression in developing panicles

We investigated gene expression changes associated with the diversity of panicle architecture and differences between the Asian and African domestication processes via RNA sequencing (RNAseq).
We used a single accession each of domesticated Asian rice (*O sativa indica* IR64) and its wild relative (*O. rufipogon* W1654), and domesticated African rice (*O. glaberrima* Tog5681) and its wild relative (*O. barthii* B88).
We also included *O. sativa japonica* cv. Nipponbare as the genomic reference accession.
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S3 **phenotype-all-varieties**).
To measure whole-transcriptome gene expression in these accessions, we collected immature panicles from each accession at four developmental stages:
rachis meristem (RM);
indeterminate meristem (IM), including panicles displaying primary branch initiation, elongation of primary branches and axillary meristem initiation;
determinate meristem (DM), including panicles where axillary meristems had transitioned into early spikelet differentiation;
and floret meristem (FM), with early differentiation of floral organs
(Figure S4a **qpcr-confirms-sampling**).
We first confirmed staging of the panicles by extracting RNA from pooled immature panicles at each stage and measuring expression of markers of panicle development by quantitative real-time RT-PCR (qPCR) (Figure S4b **qpcr-confirms-sampling**).
Because branching complexity is related to branch meristem establishment and meristem fate transition [**ref**], we chose the IM and DM stages for RNAseq.
cDNA libraries for sequencing were constucted from rRNA-depleted RNA samples from three biological replicates at both stages for all five accessions.
After mapping the reads against the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007], we obtained an average of more than 20 million uniquely mapped reads within exons for each accession, including African rice species (Table S5 **mapping-statistics**).
Pairwise distances between samples, calculated on the number of reads per gene from all detected genes, grouped samples first by stage, then by accession, and then by continent (Figure S5 **distance-heatmap**).
We did not observe grouping by domestication status, suggesting that transcriptome-wide changes during domestication are subtle compared to differences between species.

### *AP2/EREBP*-like transcription factors are core regulators of panicle branching

We used PCA on transformed raw counts to investigate general patterns of variation in the transcriptomes that may be associated with differences in panicle architecture. The first four PCs split different combinations of rice species (**supp figure?**).
PC1–PC4 may relate to species-specific differences unrelated to panicle architecture, or mapping biases introduced by mapping all libraries against the *O. sativa japonica* reference.
In contrast, PC5 separates developmental stages across all five accessions, although separation is weaker in *O. sativa indica* (Figure 3 **HB-AP2-heatmap**).
*APETALA2* and ethylene-responsive element binding protein (*AP2/EREBP*)-like genes and *MADS* genes are enriched at the extremes of PC5 (*p*~adj~ = 0.004 for both, GSEA permutation *t*-test; Figure 2 **transcriptome-pca**; Figure 3 **HB-AP2-heatmap**).
Generally, *AP2/EREBP*-like genes that contribute to PC5 are more highly expressed at the IM stage, and the *MADS* genes are more highly expressed at the DM stage (Figure 3 **HB-AP2-heatmap**).
*NAM*, *ATAF1* and *CUC2*-domain (NAC) TF genes are also more highly expressed at the IM stage, and *SQUAMOSA* promoter binding protein (SBP) genes including *WFP* and *SPL7* [**refs**] are more highly expressed at the DM stage.
There is a group of homeobox genes that is highly expressed in the DM samples, and a larger group that is highly expressed at the IM stage (Figure S6 **NAC-HB-SPL-heatmap**).

To identify the core set of genes involved in axillary meristem determination in all five accessions, we used differential expression (DE) tests to find genes that were up- or down-regulated between stages across all accessions.
Positive log~2~-fold change values (L~2~FCs) indicate higher expression in DM than in IM.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1, including 115 genes up-regulated in DM samples and 78 genes down-regulated in DM samples (Table S6 **DE-genes-stages**).
There was an enrichment of transcription factor (TF) genes in the list of 193 DE genes (38 TF genes; *p* = 1.9 × 10^-9^, hypergeometric test), including *FZP*, *LHS1*, *LAX1*, *PAP2*, and *MFO1*, which control inflorescence architecture in rice [@komatsuLAX1FRIZZYPANICLE2001; @baiRegulatoryRoleFZP2016; @baiDuplicationUpstreamSilencer2017; @kobayashiPANICLEPHYTOMER2PAP22010; @kobayashiInflorescenceMeristemIdentity2012; @khandayRiceLHS1OsMADS12013; @ohmoriMOSAICFLORALORGANS12009a], [**check Huang et al., 2018, 10.1038/ncomms3200**].
This indicates that RNAseq of developing panicles at the IM and DM stage identifies genes that control branching.
The list of DE genes included 9 *MADS* genes and 10 *AP2/EREBP*-like genes (*p*~adj~ 5.2 × 10^-11^ and 1.4 × 10^-8^ respectively, hypergeometric test; Table S6 **DE-genes-stages**).

The DE and PCA results are consistent with the role of transcriptional regulation in panicle branching, and highlight a set of candidate core regulators of axillary meristem determination and branching that are conserved in rice.
The pattern of expression of *AP2/EREBP*-like genes may indicate a role in promotion of indeterminate axillary meristem identity, or suppression of transition from axillary meristem to spikelet meristem.
*MADS* box and HB genes may have an inverse role as promoters of determinate meristem.

To test the role of *AP2/EREBP*-like genes in the control of panicle architecture, we phenotyped panicles from two loss-of-function mutants.
A mutant of the *AP2*-like gene *PLETHORA 8* [*PLT8*; @kitomiAuxinResponsiveAP22011] produces panicles with a shorter rachis with fewer primary branches (Figure 4 **panicle-mutants**), consistent with a peak of *PLT8* expression in rachis meristem tissues from *O. sativa japonica* [Figure S7 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
Panicles from *ERF142* mutants [@ayaNovelAP2TypeTranscription2014] have a reduced number of primary and secondary branches, and fewer spikelets.
*ERF142* expression is highest in primary branch and elongating primary branch meristem tissues in *O. sativa japonica* [Figure S7 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
Although neither gene was differentially expressed in our RNAseq dataset, these phenotypes support the involvement of *AP2/EREBP*-like genes in control of panicle architecture.

### *AP2/EREBP*-like gene expression is associated with panicle diversity and domestication

To identify transcriptional regulators related to variation in panicle phenotypes, we used soft clustering of scaled L~2~FCs between IM and DM to find common patterns of expression of the subset of annotated TF genes that were detected in our RNAseq dataset.
We did not include the reference accession, *O. sativa japonica* cv. Nipponbare, because it is not a high-yielding cultivar compared to other Asian rice accessions (Figure S3 **phenotype-all-varieties**).
We recovered six clusters, three of which correlated with the main principal component (PC1) in the phenotyping data (Figure 5 **cluster-phenotype-corr**).
Clusters 3, 4 and 6 correlate with spikelet number and secondary branch number but not PC1, meaning that the expression of genes in those clusters does not correlate with the phenotypic differences between wild and domesticated accessions.
Cluster 4 had an enrichment of HB genes (9 genes; *p*~adj~ = 2.5 × 10^-4^).
It also contained three MIKC^C^-type *MADS* box genes (*LHS1*, *MFO1* and *MADS14*) and three *AP2/EREBP*-like genes including *INDETERMINATE SPIKELET 1* (*IDS1*).
L~2~FCs of genes in this cluster are low in *O. sativa indica*, high in *O. glaberrima* and intermediate in the two wild accessions, suggesting divergent changes in expression pattern in the two domestications.
The correlation of L~2~FC patterns in clusters 1, 2 and 5 with PC1 suggests that genes in these clusters may be associated with changes in panicle architecture between wild and domesticated species.
Clusters 1 and 5 are positively correlated with spikelet number and secondary branch number, whereas cluster 2 is negatively correlated. 
Genes in clusters 1 and 2 had similar patterns of L~2~FC between domesticated and wild accessions in Asian and African species.
L~2~FCs are higher in domesticated accessions for genes in cluster 1, meaning that they are more highly expressed at the DM stage in domesticated accessions.
Genes in cluster 2 have higher L~2~FCs in wild accessions, meaning that they are more highly expressed at the DM stage in wild accessions.
Cluster 2 also had the strongest negative correlation with PC1, and low L~2~FCs in *O. sativa indica*. 
The lower L~2~FCs in domesticated accessions could implicate cluster 2 genes in promotion of determinate meristem fate, meaning that their lower expression at the DM stage results in more activity of indeterminate axillary meristems.
Eight *GROWTH-REGULATING FACTOR1* (GRF) family genes were detected in our dataset, and three of them were present in cluster 2 (*p*~adj~ = 7.6 × 10^-3^, hypergeometric test).
Cluster 1 and cluster 5 both had a positive correlation with PC1.
Cluster 5, which had the highest core L~2~FC in *O. sativa indica*, had an enrichment of *AP2/EREBP*-like genes (6 genes; *p*~adj~ = 7.6 × 10^-3^, hypergeometric test), and cluster 1 also contains 3 *AP2/EREBP*-like genes (Table S8 **clustered-genes**).
We used qPCR to confirm these patterns in all four stages of each accession for all *AP2/EREBP*-like genes in cluster 5 (Figure S8 **fluidigm-ap2-hb**).
Most of the genes in cluster 5 have L~2~FCs close to zero in *O. sativa indica*, and negative L~2~FCs in the other accessions (Figure S9 **cluster-5-details**), meaning that the expression of these genes decreases at the DM stage in all accessions except *O. sativa indica*.
In contrast to cluster 2 genes, the lack of repression of cluster 5 genes and to a lesser extent cluster 1 genes at the DM stage in *O. sativa indica* could result in more branching via the promotion of indeterminate/axillary meristem identity.
The enrichment of *AP2/EREBP*-like genes in cluster 5 and the presence of 3 *AP2/EREBP*-like genes in cluster 1 suggests that their pattern of expression is associated with differences in panicle architecture across wild and domesticated accessions.

To find TF genes associated with parallel changes in panicle architecture during domestication, we tested the stage × domestication interaction for *O. rufipogon*, *O. sativa indica*, *O. barthii* and *O. glaberrima* at an FDR of 0.1 (Table SX **DE-genes-interaction**).
We detected 19 genes with a stage × domestication interaction, including nine *AP2/EREBP*-like genes (*p* = 4.4 × 10^-7^, hypergeometric test; Figure 6A **dom-genes-plot**).
These genes are putative targets of parallel selection on panicle architecture that occurred during domestication.
*AP2/EREBP*-like genes were also prominent when we tested the stage × accession interaction separately for each domestication (12 out of 85 genes in Asian accessions; 8 out of 50 genes in African accessions; Table S9 **DE-genes-interaction**).
The *AP2/EREBP*-like gene *INDETERMINATE SPIKELET 1* (*IDS1*), which controls inflorescence architecture [@chuckFloralMeristemInitiation2008; @leeTwoAP2Family2012], (**check** @chengINDETERMINATESPIKELET1Recruits2018), was DE in both Asian and African domestications, although the direction of change was different (Figure 6B **dom-genes-plot**).
*IDS1* was also present in cluster 4.
Expression of genes following this pattern may have also been a target of selection on panicle architecture, but evolved divergently.

The prominence of *AP2/EREBP*-like genes among putative core regulators of branching in all four *Oryza* species, and among genes associated with differences between wild and cultivated accessions, suggest that they were key targets of artificial selection for improvement in panicle architecture, and had a role in the **reshaping of the rice transcriptome** during domestication.
