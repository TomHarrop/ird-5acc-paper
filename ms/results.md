## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice. We selected five accessions for whole-transcriptome sequencing (RNAseq). Our analysis reveals an association between expression of *AP2/EREBP*-like genes, homeobox (HB) genes and panicle architecture. We phenotyped panicles from *AP2/EREBP*-like mutants to confirm the role of those genes in branching.

### Parallel changes in panicle architecture in domesticated accessions

To measure the diversity of panicle architecture, we phenotyped 91 accessions of wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*), using `P-TRAP` image analysis software for automated measurement of traits [Figure 1A **phenotype-pca**; Figure S1 **PanicleStructure**; Table S3 **PanicleTraitsPhenotyping**; @al-tamPTRAPPanicleTrait2013].
The first principal component (PC1) in the phenotyping data accounts for 46.5% of variability.
The sum of contributions of the top four components accounted for more than 87% of the total variance, and components other than PC1 do not separate panicles from different accessions (Figure S3 **phenotype-pca-all-pc**).
PC1 separates domesticated and wild accessions, but not Asian and African accessions.
PC1 is the only component that splits the accessions by domestication status, and spikelet number, secondary branch number and primary branch number have the highest loadings on PC1 (Figure 1B **phenotype-pca**).
For all species, spikelet number correlates more with secondary branch number than it does with primary branch number.
Primary branch number correlates with spikelet number more in wild than domesticated species, but this correlation is weaker in Asian species than in African species.
Primary and secondary branch number do not correlate (**what's the cutoff for this?**), suggesting they are controlled by different genetic mechanisms (Figure S2 **correlation-pbn-spn**).
The phenotypic analysis indicates similar changes in panicle architecture between wild and domesticated accessions in the independent African and Asian domestication processes.
Spikelet number, secondary branch number and primary branch number are the main contributors to differences in panicle architecture, and these phenotypes are all related to axillary meristem formation and fate transition.
These differences may be the result of parallel, artificial selection on panicle architecture (**in discussion**?).

> TH: I'm wondering if we should add the correlation plot to figure 1, since we are relying on it a lot in this paragraph?

### Measurement of gene expression in developing panicles

We investigated gene expression differences underlying diversity of panicle architecture and differences related to the Asian and African domestication processes via RNA sequencing (RNAseq).
We used a single accession each of domesticated Asian rice (*O sativa indica* IR64) and its wild relative (*O. rufipogon* W1654), and domesticated African rice (*O. glaberrima* Tog5681) and its wild relative (*O. barthii* B88).
We also included *O. sativa japonica* cv. Nipponbare as the genomic reference accession for Asian rice.
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S4 **phenotype-all-varieties**).
To confirm the phenotypes of the accessions used for this analysis, we repeated panicle phenotyping for these accessions in **different/consistent** conditions (**to discuss: why did we do this? to control for different environments? to make sure any differences in architecture were genetic?**).
The three domesticated accessions produce more spikelets and secondary branches than their wild relatives (Figure S5 **phenotyping-mpl**; Table S4 **PanicleTraitsPhenotypingPlantsSequenced**).
The domesticated accessions have a similar number of primary branches, but Asian domesticated species have more secondary branches and spikelets than domesticated African accessions.

To measure whole-transcriptome gene expression in these accessions, we collected immature panicles from each accession at four developmental stages:
rachis meristem (RM);
indeterminate meristem (IM), including panicles displaying primary branch initiation, elongation of primary branches and axillary meristem initiation;
determinate meristem (DM), including panicles where **some/most** axillary meristems had transitioned into early spikelet differentiation;
and floret meristem (FM), with early differentiation of floral organs
(Figure S6a **qpcr-confirms-sampling**).

> I found the description of the stages complicated and confusing. I tried to make it simpler, but I'm not sure if some of the details I removed were essential. We can put more info in the figure/methods. Also, we never refer to the stages by number, so I'm not sure if we should give numbers here.

We first confirmed staging of the panicles by extracting RNA from pooled immature panicles at each stage and measuring expression of markers of panicle development by quantitative real-time RT-PCR (qPCR) (Figure S6b **qpcr-confirms-sampling**).

> Actually I don't think we even need to name the genes here, it's not the main point of the paper. We can compare them to their expected expression in the figure legend.  *LHS1*, *G1L5* (*TAWAWA*), *FZP*, *LAX1* and *MADS14*

Because branching complexity is related to branch meristem establishment and meristem fate transition [**ref**], we chose the IM and DM stages for RNAseq.
cDNA libraries for sequencing were constucted from rRNA-depleted RNA samples from three biological replicates at both stages for all five accessions.
We mapped reads against the MSU Release 7.0 of the annotation of the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007].
We observed an equivalent mapping for all the species

> **this is not really true :(** I was trying to find a tactical way of saying "we had enough mapped reads for quantification" (next sentence).

We obtained an average of more than 20 million uniquely mapped reads within exons for each accession, including African rice species (Table S5 **mapping-statistics**).
Pairwise distances between samples, calculated on the number of reads per gene from all detected genes, grouped samples first by stage, then by accession, and then by continent (Figure S7 **distance-heatmap**).
We did not observe grouping by domestication status, suggesting that transcriptome-wide changes during domestication are subtle compared to differences between species.

### Core genes related to meristem determination and/or branching


To identify the core set of genes that is involved in axillary meristem determination in all five accessions, we used differential expression (DE) tests to find genes that were up- or down-regulated between stages across all accessions.
Genes with positive L~2~FC values have higher expression in DM than in IM.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1 (Table S6 **DE-genes-stages**).
There was an enrichment of transcription factor (TF) genes in the list of 193 DE genes (38 TF genes; *p* = 1.9 × 10^-9^, hypergeometric test), including *FZP*, *LHS1*, *LAX1*, *PAP2*, and *MFO1*, which control inflorescence architecture in rice [@komatsuLAX1FRIZZYPANICLE2001; @baiRegulatoryRoleFZP2016; @baiDuplicationUpstreamSilencer2017; @kobayashiPANICLEPHYTOMER2PAP22010; @kobayashiInflorescenceMeristemIdentity2012; @khandayRiceLHS1OsMADS12013; @ohmoriMOSAICFLORALORGANS12009a], [**check Huang et al., 2018, 10.1038/ncomms3200**].
This indicates that RNAseq of developing panicles at the IM and DM stage identifies genes that control branching.
The list of DE genes included 9 *MADS* genes and 10 *APETALA2* and ethylene-responsive element binding protein (*AP2/EREBP*)-like genes (*p*~adj~ 5.2 × 10^-11^ and 1.4 × 10^-8^ respectively, hypergeometric test; Table S6 **DE-genes-stages**).
These genes, and the other DE genes, are candidate core regulators of axillary meristem determination and branching in rice.

> TH: **To discuss**. HB genes are not enriched here. Because of this and the weirdness of cluster 4, I'm worried about the emphasis on HB.

PCA on transformed raw counts also indicated that *AP2/EREBP*-like genes and homeobox (HB) genes are enriched at the extremes of a component that separates developmental stages across all five accessions (*p*~adj~ = 0.004 for both, GSEA permutation *t*-test; Figure 2 **transcriptome-pca**; Figure 3 **HB-AP2-heatmap**).
Generally, *AP2/EREBP*-like genes that contribute PC5 are more highly expressed at the IM stage (Figure 3 **HB-AP2-heatmap**).


We used PCA on transformed raw counts to investigate the patterns of variation in the transcriptomes associated with differences in panicle architecture. 
The first four PCs split different combinations of rice species (Figure 2 **transcriptome-pca**).
PC1–PC4 may relate to species-specific differences unrelated to panicle architecture, or mapping biases introduced by mapping all libraries against the *O. sativa japonica* reference.
In contrast, PC5 separates developmental stages across all five accessions, although separation is weaker in *O. sativa indica* (Figure 2 **transcriptome-pca**).
MADS, NAC and SBP genes are also enriched in PC5 (Figure S8 **NAC-MADS-SPL-heatmap**; Table S7 **PC5-TF-enrichment-gsea**).

> SJ: Indicate the objective of this to begin the paragraph. What is the additive information compared to DE analysis..(for the reader)

> TH: IMO, PCA doesn't add anything different (at the gene level) to DE analysis. It's an exploration of the results. I think we could show the component we care about (PC5) in figure 3, and move figure 2 to SI.
> **To discuss**. I think this section may be too detailed and repetitive. I've rewritten it to focus on the HB/AP2 heatmap only. 

The DE and PCA results are consistent with the role of transcriptional regulation in panicle branching.


### *AP2/EREBP*-like [and homeobox] transcription factors are differentially expressed between stages, and associated with domesticated accessions

Because of the prominence of TFs, including *AP2/EREBP*-like and HB genes, in PC5 and in DE genes between stages (Figure 3 **HB-AP2-heatmap**; Table S6 **DE-genes-stages**), we used soft clustering of scaled log~2~-fold change values (L~2~FCs) between BM and SM to find common patterns of expression of the subset of annotated TF genes that were detected in our RNAseq dataset.
We recovered six clusters.
Three clusters had low L~2~FCs in *O. sativa indica* and a negative correlation with primary branch number, secondary branch number and spikelet number.
Cluster 4 had the lowest L~2~FC in *O. sativa indica* and the strongest negative correlation with secondary branch number, and an enrichment of HB genes (9 genes; *p*~adj~ = 2.17 × 10^-5^).
Cluster 4 also contained three *AP2/EREBP*-like genes.
The genes in these clusters may promote [**SM activity**/**determinate meristem state**], and their lack of induction in *O. sativa indica* SM samples could result in more branching activity.

The other three clusters contained genes with the highest L~2~FC in *O. sativa indica* and had a positive correlation with secondary branch number and spikelet number (Figure 4 **cluster-phenotype-corr**).
These genes may promote [**axillary meristem activity**/**indeterminate meristem state**]. 
Cluster 5, which had the highest core L~2~FC in *O. sativa indica* and the highest correlations with secondary branch number and spikelet number, had an enrichment of *AP2/EREBP*-like genes (6 genes; *p*~adj~ = 7.6 × 10^-3^, hypergeometric test; Table S8 **clustered-genes**).
Most of the genes in cluster 5 have L~2~FCs close to zero in *O. sativa indica*, and negative L~2~FCs in the other accessions (Figure S9 **cluster-5-details**), meaning that the expression of these genes decreases at the SM stage in all accessions except *O. sativa indica*.
To confirm these patterns, we used qPCR on all four stages of each accession for all *AP2/EREBP*-like genes in clusters 4 and 5, and all HB genes in cluster 4 (Figure S10 **fluidigm-ap2-hb**).
**Interpret fluidigm results**.
Lack of repression in *O. sativa indica* and correlation with spikelets and higher-order branch number could mean that genes in cluster 5 specify SM, and their delayed repression in *O. sativa indica* results in more branching in this species.
The enrichment of *AP2/EREBP*-like genes in cluster 5 suggests that they are involved in this process.

To find TF genes associated with changes in panicle architecture during domestication, we tested the stage × accession interaction for African and Asian accessions separately at an FDR of 0.1 (Table S9 **DE-genes-interaction**).
For Asian accessions (*O. rufipogon* and *O. sativa indica*), there was a significant interaction for 85 genes, including 12 *AP2/EREBP*-like genes.
In African accessions, the stage × accession interaction was significant for 50 genes, including 8 *AP2/EREBP*-like genes (**check hb + ap2 numbers**).
The genes in these lists are candidate targets of artificial selection for changes in panicle architecture.
*INDETERMINATE SPIKELET 1* (*IDS1*), which controls inflorescence architecture [@chuckFloralMeristemInitiation2008; @leeTwoAP2Family2012], (**check** @chengINDETERMINATESPIKELET1Recruits2018), was DE in both comparisons.
Expression of the 10 genes whose expression is dependent on the stage × accession interaction in both species, including *IDS1* and the *AP2/EREBP*-like gene *ERF74*, may have a parallel role in the separate domestication of Asian and African rice.
**more about HB genes**.

### *AP2/EREBP* mutants have defects in panicle branching

Next, we wanted to test whether loss-of-function mutants had phenotypes consistent with our hypothesis that changes in expression of *AP2/EREBP*-like genes control panicle architecture.
A mutant of the *AP2*-like gene *PLETHORA 8* [*PLT8*; @kitomiAuxinResponsiveAP22011] produces panicles with a shorter rachis with fewer primary branches (Figure 5 **panicle-mutants**).
In laser microdissected (LMD) meristem tissues from *O. sativa japonica* panicles, *PLT8* expression peaks in the rachis meristem stage [Figure S11 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
Panicles from *ERF142* mutants [@ayaNovelAP2TypeTranscription2014] have a reduced number of primary and secondary branches, and fewer spikelets.
*ERF142* expression is highest at the primary branch and elongating primary branch stages in LMD meristem tissues [Figure S11 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
These phenotypes indicate the *AP2/EREBP*-like genes are involved in panicle architecture.
However, neither gene was differentially expressed in our RNAseq dataset, suggesting that they do not resolve the parallel evolution of panicle architecture between African and Asian rice domestications.
