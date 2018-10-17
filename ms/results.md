## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice. We selected five accessions for whole-transcriptome sequencing (RNASeq), and our analysis reveals an association between expression of *AP2/EREBP*-like genes, homeobox (HB) genes and panicle architecture. We phenotyped previously reported single-gene *AP2/EREBP*-like mutants to confirm the role of those genes in branching.

### Parallel changes in panicle architecture in domesticated accessions

To measure the diversity of panicle architecture, we phenotyped 93 accessions of wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*), using `P-TRAP` for automated measurement of traits [Figure 1A **phenotype-pca**; Table S3 **PanicleTraitsPhenotyping**; @al-tamPTRAPPanicleTrait2013].
There was a stronger correlation between secondary branch number and spikelet number in the domesticated accessions than in the wild accessions (Figure S2 **correlation-pbn-spn**).
The first principal component (PC1) in the phenotyping data accounts for 46.5% of variability.
PC1 separates domesticated and wild accessions, but not Asian and African accessions (Figure 1E **phenotype-pca**).
The sum of contributions of the top four components accounted for more than 87% of the total variance, but components other than PC1 do not separate panicles from different accessions (Figure S3 **phenotype-pca-all-pc**), and PC1 is the only component that splits the accessions by domestication status.
Spikelet number, secondary branch number and primary branch number have the highest loadings on PC1.
This indicates that these phenotypes, which are all related to axillary meristem formation and fate transition, are the main differences between panicles from wild and domesticated accessions.

### Measurement of gene expression in developing panicles

To investigate gene expression differences underlying variation in panicle architecture, we chose a single accession for each of *O. rufipogon* (W1654), *O. sativa japonica* (Nipponbare), *O. sativa indica* (IR64), *O. barthii* (B88) and *O. glaberrima* (Tog5681) for RNA sequencing (RNAseq).
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S4 **phenotype-all-varieties**), except for *O. sativa japonica* cv. Nipponbare, which we included as the reference accession for *O. sativa japonica*.
The three domesticated accessions produce more spikelets than their wild relatives (Figure S5 **phenotyping-mpl**; Table S4 **PanicleTraitsPhenotypingPlantsSequenced**).

To measure whole-transcriptome gene expression in these accessions, we collected immature panicles from each accession at four developmental stages: rachis meristem (RM); branch meristem (BM), composed of immature panicle displaying primary branch initiation, elongation of primary branches and axillary meristem initiation; spikelet meristem (SM), which includes panicles with spikelet meristem; and floret meristem (FM) (Figure S6a **qpcr-confirms-sampling**).
We first confirmed staging of the panicles by extracting RNA from pooled meristems at each stage and measuring expression of *LHS1*, *G1L5* (*TAWAWA*), *FZP*, *LAX1* and *MADS14* with quantitative RT-PCR (qRT-PCR) (Figure S6b **qpcr-confirms-sampling**).
Because branching complexity is related to branch meristem establishment and meristem fate transition [**ref**], we chose the BM and SM stages for RNAseq.
cDNA libraries for sequencing were constucted from rRNA-depleted RNA from these stages for all five accessions.
We mapped reads against the MSU Release 7.0 of the annotation of the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007].
We obtained an average of more than 20 million uniquely mapped reads within exons for each accession, including African rice species (Table S5 **mapping-statistics**).
Pairwise distances between libraries were lower for samples from the same stage, accession and continent, in that order (Figure S7 **distance-heatmap**), indicating that transcriptomic changes during domestication are subtle compared to differences between species.
We used differential expression (DE) analysis to identify genes that were up- or down-regulated between stages across all accessions.
Genes with positive L~2~FC values have higher expression in SM than in BM.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1 (Table S6 **DE-genes-stages**).
There was an enrichment of TF genes in the list of 193 DE genes (38 genes; hypergeometric test, *p* = 1.9 × 10^-9^).
DE TF genes include *FZP*, *LHS1*, *LAX1*, *PAP2*, and *MFO1*, which control inflorescence architecture in rice [@komatsuLAX1FRIZZYPANICLE2001; @baiRegulatoryRoleFZP2016; @baiDuplicationUpstreamSilencer2017; @kobayashiPANICLEPHYTOMER2PAP22010; @kobayashiInflorescenceMeristemIdentity2012; @khandayRiceLHS1OsMADS12013; @ohmoriMOSAICFLORALORGANS12009a], [**check Huang et al., 2018, 10.1038/ncomms3200**], indicating that RNAseq of developing panicles at the BM and SM stage identifies genes that control branching.

We used PCA on transformed raw counts to investigate the patterns of variation in the transcriptomes associated with differences in panicle architecture. 
The first four PCs split different combinations of rice species (Figure 2 **transcriptome-pca**).
PC1–PC4 may relate to species-specific differences unrelated to panicle architecture, or mapping biases introduced by mapping all libraries against the *O. sativa japonica* reference.
In contrast, PC5 separates developmental stages across all five accessions, although separation is weaker in *O. sativa indica* (Figure 2 **transcriptome-pca**).
The family of *APETALA2* and ethylene-responsive element binding protein (*AP2/EREBP*)-like genes are enriched at the extremes of PC5, suggesting that they contribute to phenotypic differences between BM and SM (Figure 3 **HB-AP2-heatmap**).
Homeobox, MADS, NAC and SBP genes are also enriched in PC5 (Figure 3 **HB-AP2-heatmap**; Figure S8 **NAC-MADS-SPL-heatmap**; Table S7 **PC5-TF-enrichment-gsea**).

### *AP2/EREBP*-like transcription factors are differentially expressed between stages, and associated with domesticated accessions

> When discussing enrichment in clusters, give number of genes and *p*~adj~.

Because of the prominence of transcription factor (TF) genes and *AP2/EREBP*-like genes in PC5 and in DE genes between stages (Figure 3 **HB-AP2-heatmap**; Table S6 **DE-genes-stages**), we used soft clustering of scaled log~2~-fold change values (L~2~FCs) between BM and SM to find common patterns of expression of the subset of annotated TF genes that were expressed in our RNAseq dataset.
We recovered six clusters.
Three clusters contained genes with the highest L~2~FC in *O. sativa indica*. 
These three clusters were all correlated with secondary branch number and spikelet number (Figure 4 **cluster-phenotype-corr**).
**HB cluster here**.
Cluster 5, which had the highest core L~2~FC in *O. sativa indica* and the highest correlations with secondary branch number and spikelet number, had an enrichment of *AP2/EREBP*-like genes (hypergeometric test, *p*~adj~ == **x**; Table S8 **clustered-genes**).
Most of the genes in cluster 5 have L~2~FCs close to zero in *O. sativa indica*, and negative L~2~FCs in the other accessions (Figure S9 **cluster-5-details**).
This suggests that the expression of these genes decreases in SM in all accessions except *O. sativa indica*.
To confirm this pattern, we used qPCR on all four stages of each accession for all *AP2/EREBP*-like genes in cluster 5 (Figure S10 **fluidigm-ap2-hb**).
**Fluidigm conclusion**.
Delayed or lacking repression in *O. sativa indica* and correlation with spikelets and higher-order branch number could mean that genes in cluster 5 specify SM, and their delayed repression in *O. sativa indica* results in more branching in this species.
The enrichment of *AP2/EREBP*-like genes in cluster 5 suggests that they are involved in this process.

To find TF genes associated with changes in panicle architecture during domestication, we tested the stage × accession interaction for African and Asian accessions separately at an FDR of 0.1 (Table S9 **DE-genes-interaction**).
For Asian accessions (*O. rufipogon* and *O. sativa indica*), there was a significant interaction for 85 genes, including 12 *AP2/EREBP*-like genes.
In African accessions, the stage × accession interaction was significant for 50 genes, including 8 *AP2/EREBP*-like genes (**check ap2 numbers**).
The genes in these lists are candidate targets of artificial selection for changes in panicle architecture.
*INDETERMINATE SPIKELET 1* (*IDS1*), which controls inflorescence architecture [@chuckFloralMeristemInitiation2008; @leeTwoAP2Family2012], (**check** @chengINDETERMINATESPIKELET1Recruits2018), was DE in both comparisons.
Expression of the 10 genes whose expression is dependent on the stage × accession interaction in both species, including *IDS1* and the *AP2/EREBP*-like gene *ERF74*, may have a parallel role in the separate domestication of Asian and African rice.

### *AP2/EREBP* mutants have defects in panicle branching

Next, we wanted to test whether loss-of-function mutants had phenotypes consistent with our hypothesis that changes in expression of *AP2/EREBP*-like genes control panicle architecture.
Mutant generation in rice takes up to **x** months/years (**reference**), so we chose to start by phenotyping panicles in mutants that were publicly available.
A **knockout?** mutant of *PLETHORA 8* (*PLT8*), an *AP2*-like gene characterised in (**ref**), produces a shorter rachis with fewer primary branches than the background accession (Figure 5 **panicle-mutants**), consistent with its reported peak in expression in rachis meristem tissue in *O. sativa japonica* cv. Nipponbare (Figure S11 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016).
*ERF142* and *DLT* mutants, reported in (**ref**), produce fewer primary and secondary branches and fewer spikelets (Figure 5).

>**in discussion**: We detected all three genes at both stages of all five accessions in our RNAseq dataset, but none of them were differentially expressed. Although this suggests that they were not targets of domestication, the phenotypes support a role for AP2/EREBP-like genes in panicle architecture.
> Could we say these genes are involved in branching but not necessarily domestication - that's why they didn't show up in our analysis?
