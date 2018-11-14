## Figure legends

**Figure 1 phenotype-pca**.
The main component of variablity in panicle phenotypes from 91 rice accessions splits accessions by domestication status, and is related to spikelet number, secondary branch number and primary branch number.
(**A**) We measured traits using spread panicles from *O. rufipogon*, *O. sativa*, *O. barthii* and *O. glaberrima*.
(**B**) The first principal component (PC1) in the panicle phenotye data accounts for 46.5% of variability and separates wild and domesticated accessions.
(**C**) Spikelet number (SpN), secondary branch number (SBN) and primary branch number (PBN) has the highest loadings on PC1.
(**D**) Correlation between the main panicle traits that contribute to panicle architecture diversity.
Primary branch number and spikelet number correlate in wild species.
Secondary branch number and spikelet number correlate more in cultivated species than in wild species.
Primary and secondary branch numbers don't correlate.
RL: Rachis length; PBL: Primary branch length; PBIL: Primary branch internode length; SBL: Secondary branch length; SBIL: Secondary branch internode length; TBN: Tertiary branch number; PL: Panicle length.

**Figure 2 transcriptome-pca**.

> Move to SI

Principal components analysis of transformed read counts for each library.
Principal component 5 (PC5) separates RNAseq samples by developmental stage, and explains 5.4% of total variability.
The first four components explain 51.7% of variability, and separate RNAseq samples by species.

>TH. Reformat: Change the labels on the x-axis to match the stages in the paper. Remove the species codes from the x-axis labels & use full species names in the facet labels. Remove the guide for the bar colour. Use brackets rather than ~ to separate PC number from percent variability.

**Figure 3 HB-AP2-heatmap**.
*AP2/EREBP*-like and homeobox (HB) transcription factors change expression between BM and SM.
**A** *AP2/EREBP*-like and HB genes are distributed at the extremes of genes ranked on PC5 (*p*~adj~ = 0.004 for both, GSEA permutation *t*-test).
For the heatmap, we used the 10% of genes that have the highest absolute loading on PC5, shown in red in the enrichment plot.
**B** Most *AP2/EREBP*-like genes that pass the cutoff are more highly expressed in the BM.
Three of the four *AP2/EREBP*-like genes that are more highly expressed in the SM belong to the *AP2* subfamily.
Genes that are more highly expressed in BM mainly belong to RAV, DREB and ERF subfamilies.
**C** Most HB genes that pass the cutoff are more highly expressed in the SM.
Ten out of twenty of these genes belong to the HD−ZIP IV subfamily.

> TH: Otho, can you explain what's on the scale bar in figure 3B? Do we need a citation for the permutation test? Do you have any stats for Figure 3C? I could do a hypergeometric test if you like.

**Figure 4 cluster-phenotype-corr**.
Gene expression clusters correlate with the main component of divesity of panicle architecture (PC1), the number of secondary branches (SBN) and spikelets (SpN).
(**A**) Mean, scaled log~2~-fold change (L~2~FC) of genes by cluster and accession.
(**B**) Pearson correlation with PC1, SBN and SpN.
PC1 is the main principal component that separates domesticated and wild accessions of rice (**Figure 1**).
Correlations with SBN and SpN are based on repeat panicle phenotyping for the accessions used for RNAseq in greenhouse conditions (Figure S5 **phenotyping-mpl**; Table S4 **PanicleTraitsPhenotypingPlantsSequenced**).]
The genes in clusters 1 and 2 may be related to domestication, because L~2~FC of genes in those clusters are higher or lower in domesticated accessions, respectively, and both clusters correlate with PC1 and SpN.
Clusters of genes with high L~2~FC in *O. sativa indica* have a positive correlation with SBN and SpN, whereas clusters of genes with low L~2~FC in *O. sativa indica* have a negative correlation with SBN and SpN.
Six of the 17 genes in Cluster 5, which has the highest L~2~FC in *O. sativa indica* and the highest correlation with SBN, are *AP2/EREBP*-like genes (*p*~adj~ = 7.6 × 10^-3^, hypergeometric test; Table S8 **clustered-genes**).
Conversely, nine of the 31 genes in Cluster 4, which has the lowest L~2~FC in *O. sativa indica* and the strongest negative correlation with SBN, are homeobox genes (*p*~adj~ = 2.17 x 10^-5^, hypergeometric test; Table S8 **clustered-genes**).

**Figure 5 panicle-mutants**.

> To move to SI?

Mutants in two *AP2/EREBP*-like genes, *PLT8* and *ERF142*, have defects in panicle architecture compared to their background accessions.
The *PLT8* mutant produces fewer primary branches and spikelets, and the *ERF142* mutant produces fewer primary branches, secondary branches and spikelets.

> TH: Otho, can you use the official gene names for the figures? *PLT8* and *ERF142*.

**Figure 6 dom-genes**.
Parallel and divergent evolution of gene expression during domestication.
(**A**) Expression of genes with a stage × domestication interaction when both domestications are considered have parallel changes in expression pattern in panicles at indeterminate (IM) and determinate (DM) meristem stage.
(**B**) Some of the genes with a stage × accession interaction in both domestications, when African and Asian accessions are tested separately, have divergent changes in expression between wild and domesticated accessions.

## Supplementary figures

**Figure S1 PanicleStructure**.
Spread mature rice panicle. RL: Rachis length; PB: Primary branch; PBL: Primary branch length; PBintL: Primary branch internode length; SB, Secondary branch; SBL: Secondary branch length; SBintL: Secondary branch internode length; Sp : Spikelet.

**Figure S2 correlation-pbn-spn**.
Moved to figure 1D.

**Figure S3 phenotype-pca-all-pc**.
Principal component analysis (PCA) of panicle phenotyping data showing components 1–4.
PC1 accounts for 46.5% of variability and separates panicles from domesticated and wild accessions.
The lower ordinates do not separate panicles by species.

**Figure S4 phenotype-all-varieties**.
The accessions used for RNAseq are consistent with species-wide patterns of panicle architecture.
The *y*-axis shows the projection of each panicle on principal component 1 (PC1), which separates wild and domesticated accessions (**Figure 1**).
We measured 3 to 9 panicles from 1 to 3 plants for each accession.
Scores on PC1 for the accessions chosen for RNAseq are shown in red.

**Figure S5 phenotyping-mpl**.
Detailed phenotyping of the five Oryza accessions used for RNAseq.
We phenotyped 9 panicles from each accession.
Plants for this dataset were grown together in controlled, greenhouse conditions, separately to the plants used for the 91-accession survey, using the same conditions as plants that were used for RNAseq.
The three domesticated accessions produce more spikelets and secondary branches than their wild relatives.
The domesticated accessions have a similar number of primary branches, but Asian domesticated species have more secondary branches and spikelets than domesticated African accessions.
In comparison to *O. sativa japonica*, *O. sativa indica* produces more secondary branches.

**Figure S6 qpcr-confirms-sampling**.
Early stages of rice panicle development used for gene expression analysis.
(**A**) Developmental stages of immature panicles collected for expression analysis.
Stage 1: rachis meristem;
Stage 2: indeterminate meristem (IM) stage with formation of primary branch meristems, elongation of primary branch meristem and formation of axillary meristem;
Stage 3: determinate meristem (DM) stage with spikelet meristem and floret differentiation;
Stage 4: floret displaying early floral organ differentiation.
The scale bar indicates 100 μm.
(**B**) Quantitative RT-PCR using meristem stage-specific marker genes for validation of staging.
*LHS1* is used as a marker of spikelet meristem differentiation [**ref**].
*FZP* is used as… LAX1….TAW …. osMADS14…. [**refs**].
AM, axillary meristem;
SpM, spikelet meristem;
RM, Rachis meristem;
PBM, primary branch meristem;
ePBM, elongating primary branch meristem;
FlM, floret meristem;
St, stamen;
p, palea;
l, lemma.

> TH: Relabel to match stage abbreviations from results.
> Change PbM to PBM and ePbM to ePBM

**Figure S7 distance-heatmap**.
Heatmap of pairwise distances between RNAseq samples. Samples group by stage, species and continent. ob, *O. barthii*; og, *O. glaberrima*; osj, *O. sativa japonica*; osi, *O. sativa indica*; or, *O. rufipogon*; PBM, branch meristem; SM, spikelet meristem

> TH: (fixme) Relabel axes to make them human-readable and to match abbreviations from results. Remove replicate numbers.

**Figure S8 NAC-MADS-SPL-heatmap**.
Expression of genes from selected transcription families.
While *MADS* and *SBP* genes are expressed preferentially in the SM, *NAC* are more expressed in the BM.
(**A**) Enrichment is estimated with the the GSEA method, which returns
a permutation based adjusted *p*-value of 0.0037 for MADS box genes, of 0.0.0055 for NAC genes and of 0.0245 for SBP genes.
(**B−D**) For the heatmap, we used the 10% of genes that have the highest absolute loading on PC5 (shown in red in the enrichment plot).
Heatmaps display normalized RNAseq counts which have been z−score scaled independently for each species.

>HA : BM for IM and SM to DM in heatmap and on the pC5 value (x-label). Species in italique ? Just add  « B » for the second panel in this figure. (remove « C »). Order of the species – Remove text under figure

**Figure S9 cluster-5-details**.
Most genes in cluster 5 have negative L~2~FCs between BM and SM in *O. rufipogon*, *O. barthii* and *O. glaberrima*, but L~2~FCs in *O. sativa indica* are closer to zero.
This cluster has an enrichment of *AP2/EREBP*-like genes.

> TH: Fixme. Relabel to match stage abbreviations from results. Add cluster 4?

**Figure S10 fluidigm-ap2-hb**.
Expression analysis along early panicle development of AP2-like and Homeobox genes present in Cluster 4 and 5.
Stage 1 to 4 corresponds respectively to Rachis Meristem, Indeterminate Meristem, Determinate meristem and Flower.

> TH: Relabel to match stage abbreviations from results, and remove japonica.
> 
> TH: I'm not sure if the qPCR helps for cluster 4. Most of the time, indica and rufipogon look the same (to me). 

**Figure S11 lmd-paper-ap2**.
Expression of *AP2/EREBP*-like genes in *O. sativa japonica* cv. Nipponbare meristems [data from @harropGeneExpressionProfiling2016].
Both genes are expressed at all stages.
*PLT8* expression peaks in RM.
RM, rachis meristem; PBM, primary branch meristem; ePBM/AM, extending primary branch meristem and axillary meristem; SM, spikelet meristem.

> TH: (fixme) make this look nicer

## Supplementary tables

**Table S1 Plantinfo**.
Detailed information of rice accessions used in this study.

**Table S2 supp-table-PrimerList**.
Sequences of primers used.

**Table S3 PanicleTraitsPhenotyping**.
Detailed quantification of panicle traits in 93 accessions from wild and domesticated Asian and African rice species.
ob, *O. barthii*;
og, *O. glaberrima*;
os, *O. sativa*;
or, *O. rufipogon*;
RL, Rachis length;
PBN,Primary branch;
PBL, Primary branch length;
PBIL, Primary branch internode length;
SBN, Secondary branch;
SBL,Secondary branch length;
SBIL, Secondary branch internode length;
SpN, Spikelet number;
TBN, Tertiary branch number.

**Table S4 PanicleTraitsPhenotypingPlantsSequenced**.
Detailed quantification of panicle traits from rice accessions used for sequencing analysis.
These plants were grown and the same time and same conditions as the plants used for gene expression analysis. 

**Table S5 mapping-statistics**.
Read and mapping statistics for all RNAseq samples.

**Table S6 DE-genes-stages**.
Differential expression test results between stages across all species.

**Table S7 PC5-TF-enrichment-gsea**.
Transcription factor families that are enriched along PC5.

**Table S8 clustered-genes**.
Clustered genes.

**Table S9 DE-genes-interaction**.
Differential expression test results for the stage × accession interaction in Asian and African accessions.

**Table S10 PanicletraitsPhenotypingAP2**.
Detailed quantification of panicle traits from crl5 and smos1 mutants.
RL, Rachis length;
PBN,Primary branch;
PBL, Primary branch length;
PBIL, Primary branch internode length;
SBN, Secondary branch;
SBL,Secondary branch length;
SBIL, Secondary branch internode length;
SpN, Spikelet number;
TBN, Tertiary branch number

> TH: gene names