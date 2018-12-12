## Figure legends

**Figure 1 phenotype-pca**.
The main component of variablity in panicle phenotypes from 91 rice accessions splits accessions by domestication status, and is related to spikelet number, secondary branch number and primary branch number.
(**A**) We measured traits using spread panicles from *O. rufipogon*, *O. sativa*, *O. barthii* and *O. glaberrima*.
(**B**) The first principal component (PC1) in the panicle phenotype data accounts for 46.5% of variability and separates wild and domesticated accessions independently to continent.
(**C**) Spikelet number (SpN), secondary branch number (SBN) and primary branch number (PBN) have the highest loadings on PC1.
(**D**) Correlation between the main panicle traits that contribute to panicle architecture diversity.
Primary branch number and spikelet number correlate in wild species.
Secondary branch number and spikelet number correlate more in cultivated species than in wild species.
Primary and secondary branch numbers do not correlate.
RL: Rachis length; PBL: Primary branch length; PBIL: Primary branch internode length; SBL: Secondary branch length; SBIL: Secondary branch internode length; TBN: Tertiary branch number.

**Figure 2 transcriptome-pca**.
Principal components analysis of transformed read counts for each library.
Principal component 5 (PC5) separates RNAseq samples by developmental stage, and explains 3.5% of total variability.
The first four components explain 86.3% of variability, and separate RNAseq samples by species.
Bars show single samples (three replicates per accession per stage).

**Figure 3 MADS-AP2-heatmap**.
*AP2/EREBP*-like and MADS-box transcription factors change expression between IM and DM.
For each gene family, we plotted the 10% of genes that have the highest absolute loading on PC5.
Genes in the upper panels have a positive loading on PC5 (corresponding to higher expression at the IM stage), whilst genes in the lower panels had a negative loading.
Most *AP2/EREBP*-like genes that pass the 10% cutoff are more highly expressed in the IM.
Two of the four *AP2/EREBP*-like genes that are more highly expressed in the DM belong to the AP2 clade.
*AP2/EREBP*-like genes that are more highly expressed in IM mainly belong to the ERF clade.
In contrast, most MADS-box genes that pass the cutoff are more highly expressed at the DM stage.
Clades for *AP2/EREBP*-like genes are from the Plant Transcription Factor Database v4.0 [@jinPlantTFDBCentralHub2017].
MADS-box clades were manually tabulated from Figure 2 of @aroraMADSboxGeneFamily2007 (**Table S11 arora_clades**).

**Figure 4 panicle-mutants**.
Mutants in two *AP2/EREBP*-like genes, *PLT8* and *ERF142*, have defects in panicle architecture compared to their background accessions.
The *crl5* mutant of *PLT8* (*LOC_Os07g03250*) produces fewer primary branches and spikelets, and the *smos1-3* mutant of *ERF142* (*LOC_Os05g32270*) produces fewer primary branches, secondary branches and spikelets.

**Figure 5 cluster-phenotype-corr**.
Gene expression clusters correlate with the main component of divesity of panicle architecture (PC1), the number of secondary branches (SBN) and spikelets (SpN).
Clusters contained 19–31 genes each (**Table S8**).
(**A**) Mean, scaled log~2~-fold change (L~2~FC) of genes by cluster and accession.
(**B**) Pearson correlation with PC1, SBN and SpN.
PC1 is the main principal component that separates domesticated and wild accessions of rice (**Figure 1**).
Correlations with SBN and SpN are based on repeat panicle phenotyping for the accessions used for RNAseq in greenhouse conditions (Figure S10 **phenotyping-mpl**; Table S4 **PanicleTraitsPhenotypingPlantsSequenced**).]
The genes in clusters 1 and 2 may be related to domestication, because L~2~FC of genes in those clusters are higher or lower in domesticated accessions, respectively, and both clusters correlate with PC1 and SpN.
Clusters of genes with high L~2~FC in *O. sativa indica* have a positive correlation with SBN and SpN, whereas clusters of genes with low L~2~FC in *O. sativa indica* have a negative correlation with SBN and SpN.
Six of the 17 genes in Cluster 5, which has the highest L~2~FC in *O. sativa indica* and the highest correlation with SBN, are *AP2/EREBP*-like genes (*p*~adj~ = 7.6 × 10^-3^, hypergeometric test; Table S8 **clustered-genes**).
Conversely, nine of the 31 genes in Cluster 4, which has the lowest L~2~FC in *O. sativa indica* and the strongest negative correlation with SBN, are homeobox genes (*p*~adj~ = 2.17 x 10^-5^, hypergeometric test; Table S8 **clustered-genes**).

**Figure 6 dom-genes**.
Parallel and divergent evolution of gene expression during domestication.
(**A**) Expression of genes with a stage × domestication interaction when both domestications were considered have parallel changes in expression pattern in panicles at indeterminate (IM) and determinate (DM) meristem stage.
(**B**) Expression of genes with a stage × accession interaction in both domestications when tested separately, but not a stage × domestication interaction when both domestications were tested together.
These genes have divergent changes in expression between wild and domesticated accessions.

## Supplementary figures

**Figure S1 PanicleStructure**.
Spread mature rice panicle. RL: Rachis length; PB: Primary branch; PBL: Primary branch length; PBIL: Primary branch internode length; SB, Secondary branch; SBL: Secondary branch length; SBIL: Secondary branch internode length; Sp : Spikelet.

**Figure S2 phenotype-pca-all-pc**.
Principal component analysis (PCA) of panicle phenotyping data showing components 1–4.
PC1 accounts for 46.5% of variability and separates panicles from domesticated and wild accessions.
The lower ordinates do not separate panicles by species.

**Figure S3 phenotype-all-varieties**.
The accessions used for RNAseq are consistent with species-wide patterns of panicle architecture.
The *y*-axis shows the projection of each panicle on principal component 1 (PC1), which separates wild and domesticated accessions (**Figure 1**).
We measured 3 to 9 panicles from 1 to 3 plants for each accession.
Scores on PC1 for the accessions chosen for RNAseq are shown in red.

**Figure S4 qpcr-confirms-sampling**.
Early stages of rice panicle development used for gene expression analysis.
(**A**) Developmental stages of immature panicles collected for expression analysis.
Stage 1: rachis meristem;
Stage 2: indeterminate meristem (IM) stage with formation of primary branch meristems, elongation of primary branch meristem and formation of axillary meristem;
Stage 3: determinate meristem (DM) stage with spikelet meristem and floret differentiation;
Stage 4: floret displaying early floral organ differentiation.
The scale bar indicates 100 μm.
(**B**) Quantitative RT-PCR using meristem stage-specific marker genes for validation of staging.
AM, axillary meristem;
SpM, spikelet meristem;
RM, Rachis meristem;
PBM, primary branch meristem;
ePBM, elongating primary branch meristem;
FlM, floret meristem;
St, stamen;
p, palea;
l, lemma.

**Figure S5 distance-heatmap**.
Heatmap of pairwise distances between RNAseq samples.
Samples group by stage, species and continent.
The numbers indicate single samples (three replicates per accession per stage).
The axes are ordered by hierarchical clustering of Minkowski distances between samples.

**Figure S6 NAC-HB-SPL-heatmap**.
Expression of genes from selected transcription families.
While Homeobox and *SBP* genes are expressed preferentially in the DM, *NAC* are more expressed in the IM.
We used the 10% of genes that have the highest absolute loading on PC5.
Genes in the upper panels have a positive loading on PC5 (corresponding to higher expression at the IM stage), whilst genes in the lower panels have a negative loading.
Clades for Homeobox genes are from the Plant Transcription Factor Database v4.0 [@jinPlantTFDBCentralHub2017].

**Figure S7 lmd-paper-ap2**.
Expression of *AP2/EREBP*-like genes in *O. sativa japonica* cv. Nipponbare meristems [data from @harropGeneExpressionProfiling2016].
Both genes are expressed at all stages.
*PLT8* expression peaks in RM.
RM, rachis meristem; PBM, primary branch meristem; ePBM/AM, extending primary branch meristem and axillary meristem; SM, spikelet meristem.

**Figure S8 fluidigm-ap2-hb**.
Expression analysis along early panicle development of AP2-like and Homeobox genes present in Cluster 4 and 5.
Stage 1 to 4 corresponds respectively to Rachis Meristem, Indeterminate Meristem, Determinate meristem and Flower.

**Figure S9 cluster-5-details**.
Most genes in cluster 5 have negative L~2~FCs between BM and SM in *O. rufipogon*, *O. barthii* and *O. glaberrima*, but L~2~FCs in *O. sativa indica* are closer to zero.
This cluster has an enrichment of *AP2/EREBP*-like genes.

**Figure S10 phenotyping-mpl**.
Detailed phenotyping of the five Oryza accessions used for RNAseq.
We phenotyped 9 panicles from each accession.
Plants for this dataset were grown together in controlled, greenhouse conditions, separately to the plants used for the 91-accession survey, using the same conditions as plants that were used for RNAseq.
The three domesticated accessions produce more spikelets and secondary branches than their wild relatives.
The domesticated accessions have a similar number of primary branches, but Asian domesticated species have more secondary branches and spikelets than domesticated African accessions.
In comparison to *O. sativa japonica*, *O. sativa indica* produces more secondary branches.

## Supplementary tables

**Table S1 Plantinfo**.
Rice accessions used in this study.

**Table S2 supp-table-PrimerList**.
Sequences of primers used.

**Table S3 PanicleTraitsPhenotyping**.
Quantification of panicle traits in 91 accessions from wild and domesticated Asian and African rice species.
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
These plants were grown at the same time and in the same conditions as the plants used for gene expression analysis. 

**Table S5 mapping-statistics**.
Read and mapping statistics for all RNAseq samples.

**Table S6 DE-genes-stages**.
Differential expression test results between stages across all species.

**Table S7 PC5-TF-enrichment-gsea**.
Transcription factor families that are enriched along PC5.

**Table S8 clustered-genes**.
Clustered genes with scaled L~2~FC values for each accession.

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

**Table S11 arora_clades**.
Clade information for MADS-box genes, manually tabulated from Figure 2 of @aroraMADSboxGeneFamily2007.
