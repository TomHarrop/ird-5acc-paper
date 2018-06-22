---
title: "Temporary - Determination of Inflorescence Architecture in Wild and Domesticated Rices"
bibliography: ../bib/references.bib
---
# Introduction

* evo devo and gene expression
* transcription factors / evolution of gene regulation
* domestication and gene expression (maize, tomato transcriptomes)
* branching development in rice : undeterminated vs determinate meristem
* details about rice domestication

* **aims**:
    - Genes/pathway involved in the inflorescence determination in rice  
    - changes in phenotype, gene expression caused by / associated changes/switch/delay in determination of inflorescence
    - Changes in phenotype link to domestication
    - common genes that respond to selection in different domestications
    - evolution of transcription factor family expression
    - does selection for the same phenotype affect expression of the same genes / sets of genes


# Results

## Panicle traits of independently domesticated rice species are similar

We have measured panicle traits of 91 African and Asian rice accessions ***(Supplementary Dataset 1)***. We show that the main traits that differ among accessions are branch number and spikelet number, and that those traits split wild and domesticated species.

Those 91 rice accessions belong to 4 rice species: the wild African rice *Oryza barthii*, the wild Asian rice *Oryza Rufipogon*, the domesticated African rice *Oryza glaberrima* and the domesticated Asian rice,  *Oryza sativa*. We have measured ***approximatively???*** 20 accessions per species and 10 plants per accession for a total of 1140 ***panicles/plants???***.

For each plant, we have measured these phenotypic traits using the software pTrap as reported by [@faroq2013p]:

- Rachis Length (RL),
- Primary Branch Numbers (PbN),
- Primary Branch Average Length (PbL),
- Primary Branch Internode Average length (SbIntl),
- Secondary Branch Numbers (SbN),
- Secondary Branch Numbers (SbIntl),
- Secondary Branch Average Length (SbL),
- Tertiary Branch Numbers (TbN),
- Spikelet Numbers (SpN).

The recorded traits are in Supplementary Table 1.

In that dataset, we have searched for recurrent patterns of phenotypic traits using Principal Component (PC) Analysis ***(Fig. 01)***. The PC analysis shows that most of the variability in the dataset is explained by the domestication factor, because the first PCs, that collects almost half of traits variability, perfectly splits wild and domesticated rice species. Interestingly, those domesticated species belong to both African and Asian rice lineages, which have been domesticated independently. Indeed, independent lineage seems to have little effects on the results of PC analysis, suggesting that the domestication has a much greater effect on panicle traits the origin lineage of different accessions.

The rotations of the first PC ***(Fig. 01)*** shows that all panicle traits have adapted have slightly adapted with domestication. Indeed those rotations split traits that are expected to increase with spikelet number (PbN, SbN, RL, PbL, SbL, PanL and TbN) from those that are expected to decrease with spikelet number (PbIntl and SbIntL).

## PbN and SbN are independent and they both correlate with SpN

We have used our trait dataset to explore the phenotype of wild and domesticated rice panicles. We concentrated on the traits that influence most the first PC: PbN, SbN and SpN ***(Fig. 2)***.

SpN is highly different between wild and domesticated species; with wild rice accessions producing between 6 and 186 spikelets per panicle and cultivated accessions producing between 62 and 391 spikelets per panicle ***(Fig. 2)***.

The trait that correlates most with SpN is SbN (Fig. 2). Indeed the linear correlation coefficients between the two, in the four species, range between 0.83 and 0.93. The lowest correlation coefficient belongs to *O. rufipogon*. Indeed, a close inspection of the scatter plots of this species reveals that at low SpN, the linear relationship between SbN and SpN is lost. Probably because many accessions of *O. rufipogon* rely more heavily on primary branching.

PbN correlates with SpN with coefficients ranging between 0.33 and 0.86. The lowest coefficient belongs to *O. sativa*, probably because this species produces the highest and the most variation in SbN. A similar behaviour is displayed by *O. glaberrima*: the other domesticated species.

Primary and Secondary branch number correlate poorly, with correlation coefficients ranging between 0.30 and 0.63. The lowest coefficients belong to domesticated species. This suggests that SpN and SbN are controlled by distinct genetic mechanism.

## We have selected 5 rice accessions for RNAseq

To investigate what determines branching, we have selected five rice accessions. Those accessions were selected because they grow easily in greenhouse and/or because they carry distinctive traits. Details on those 5 species are in ***table 1***.  

> **Details on RNAseq here? Sampling, number of reads, mapping efficiency**

On the same exact plants that we have used for RNAseq, we have measured panicle phenotypic traits ***(Supplementary Dataset 2)***.  Those phenotypic traits match what expected from earlier measurements ***(Fig-03)***.

## Transition to the Spikelet Meristem is sharper in low branching species

We have compared the transcriptomes of PBM and SM separately in those 5 rice species and we have detected with a significance level of alpha = .05:

- 1942 differentially expressed genes in ***O. barthii***.
- 1848 differentially expressed genes in ***O. glaberrima***.
- 742 differentially expressed genes in ***O. rufipogon***.
- 387 differentially expressed genes in ***O. indica***.
- 831 differentially expressed genes in ***O. japonica***

> **The number of DE genes is lower in species that make more SbN. How do we discuss this?**

## Transcriptomic patterns in multispecies experiments can be detected with PCA

In RNAseq experiments that compare (multiple conditions from) multiple species, useful and reliable expression pattern are often layered under various sources of noise or of unreliable information: that noise is a mix of mapping bias and species specific expression values that are unrelated to the conditions of interest.

We have successfully used principal component analysis (PCA) to detect and remove those layers of unreliable information and to reveal expression patterns that are concealed underneath them.

We have achieved this by running PCA on rlog transformed and normalized counts [***rlog deseq2 citation***] ***(Fig-04)***.

The rotations of first 4 PC split different combinations of rice species: PC1 splits African and Asian species, PC2 splits *O. Nipponbare* (which was used as mapping reference) from the two other Asian species, PC3 splits *O. indica* from *O. rufipogon*, PC4 splits *O glaberrima* from *O. barthii*. These four component do not enrich functional categories. Possibly, because they mix differential expression with mapping bias and unrelated species specific differences.

The PC5 instead splits developmental stages across all species but *O. indica*, and also enriches functional categories. We have decided to explore the PC5 more in detail

## Multiple transcription factor families are co-expressed in the branch meristems across species

We have used principal component (PC) analysis to group genes that behave similarly in the various species. The fifth PC (PC5) successfully groups genes that are co-expressed across the species. Many transcription factor families are enriched among those genes ***(Fig-05, Suppl-Fig-01)***.

### AP2

AP2 genes [***here a good review of ap2s?***] are the most differentially expressed between BM and SM. While only 2 AP2 genes are upregulated in the SM, at least 20 AP2 genes are upregulated in the BM.

> ***Extend***

### MADS

Eight MADS-box genes, many of which are confirmed to be involved in spikelet development, are upregulated in the SM. Among them *AFO*, *PAP2*, *OsMADS15 - DEP*, *OsMADS55*, *OsMADS6 - MFO*, *OsMDS14*, *OsMADS17*, *OsMADS32 - CFO1*. Other four could be over-expressed specifically in African species: *OsMADS4*,*OsMADS45*, *OsMADS20*, *OsMADS2*, *OsMADS18*.

Five MADS-box genes are upregulated instead in the BM. Among them: *MADS16 - SPW1*, *LOC_Os06g01890*, *OsMADS47 - OsMDP1* [root elongation], *OsMADS26*, *OsMADS50*

> ***Reshape, extend***

### zf-HD

Six *zf-HD* genes are upregulated in the SM, especially in African species. Among them *osZHD4* and *osZHF8*.
One zf-HD gene is upregulated in the BM, also especially in African species.

zf-HD proteins seems to have a specific function in African Species.

> ***Extend***

### Homeobox

Several *Homeobox* are upregulated in the SM. Among them *LSY1* (regulator of adaxial/abaxial? trichome? interacts with AP2s) and several *HB - START domain* (They interact with Lipids - also several Lipid Binding Proteins are Upregulated, or maybe same annotations for the same genes).

Also one *Homeobox* gene is upregulated in the BM of African species: *LOC_Os10g26500*

> ***Extend, revise: get better functional annotation for HB subfamilies***

### SPB

At least four *SPB* genes are upregulated in the SM: *SPL10*, *SPL13*, *SPL8-SGL1*, *SPL7*, *SPL10*.

*SPL 3*, instead is upregulated in the BM.

### NAC

At least eight *NAC* genes are upregulated in the SM: among the *NAC9*, *NAC3*, *NAC14*, *NAC22*, *OMTN1*; (some are involved in root system architecture).

. **Presentation of the stage/time course, checking some genes by Q-PCR (perhaps other genes necessary which could be include in the next fluidgm)**

> **What are the genes/patwhay involved in the meristem state in rice inflorescence Undeterminate vs determinate**



### Gene families related to differentiation are enriched in the SM

Those genes that are expressed preferentially in the SM, rank negatively in the loadings 5^th^ PC. Those genes are often linked to organ differentiation.

Among those, 8 MADS-box genes are among the 100 genes that rank most negatively.

![](../fig/fig-tmp-MADS-bot_pc5.svg)

- 7 of them: *OsMADS32* [@sang2012chimeric], *OsMADS6* [@ohmori2009mosaic] *OsMADS1* [@prasad2005osmads1], *OsMADS14* and *OsMADS15* [@wu2017abcs], *OsMADS34* [@kobayashi2009panicle], *OsMADS5* [@cui2010functional] determine floral organ identity and spikelet development.

- 1 of them, *LOC_Os04g49150* is uncharacterized and it is a new candidate for the same role as the genes mentioned above.

> **Although it is expressed at lower level? can it be a mapping artifact?**

Also three YABBY genes are expressed preferentially in the SM: *OsYABBY4* [@yang2016rice], *OsYABBY5* and *OsYABBY3*. Among those *OsYABBY5* determines the identity of the spikelet meristem [@tanaka2012yabby].

### PBM - Top PC5 - Indica same

![](../fig/fig-TMP-pc5-top200-cl3.jpeg)

Interesting genes highlighed : AP2/ALOG on one side and MADS box genes in the other side of the FigPCA. There is also genes related to hormone.

Interestingly, some genes seems to have a a variant pattern in indica (some AP2 and MADS box) : these genes are present on cluster3 and 4
Maybe make a bridge with phenotype (high and low branching) that could be explained by these genes.

### Cluster analysis : This part have to be update with the last clustering data

> **the different pattern could be link to phenotype (correlation test) and or are interesting based on the gene enrichment and composition**

![fig-cluster: TF cluster ](../fig/cluster.pdf)
Cluster3 : high branching phenotype; divergence of slope for indica (positive slope in indica, negative for the other). This cluster is enriched in AP2 genes (gsea analysis) and most of these genes are positively correlated with Sb/sp nb
Cluster 4 : branching phenotype; divergence of slope for indica (negative slope in indica), this cluster is enriched in HB TF, are present also some genes known to be involved in sp determination (LHS1, MADS14, MFO1, OsIDS1). Enriched also in SWI/SNF-SW
The HB are not negatively correlated with SP, but it seems to have a pattern for other genes.

Cluster 1 : related to domestication/Phenotype; positive slope in wild species vs negaive slope in domesticated
This cluster is enriched in genes related to chromatin regulation (PHD, CCP, SNF2..).
These genes are negatively correlated with SpNb.
Cluster6 : inverse pattern to cluster 1 and positive correlation with SpNB...At this time, i don't see clear interesting genes  enrichisment

Cluster 5 : related to domestication. For both domestication process, we observe a change of slope between wild and domesticated but with different level between asian and african.
In this cluster : LAX1, APO2, MOC1, SHAT5
no correlation with Phenotype
Link with PC3 and PC4 which split the wild/domesticated in asian and african separatly

Cluster 2 : related to continent?
different slope between asian and african... find a link with PCA analysis and/or DE analysis (stage~continent)


The difficulties is to represent the merging of the data all together : how to do it?

> **Next work to finalize
Analyse along time course AP2/HB genes presents in cluster 3 and 4 to evaluate if there is a delay in expression that could explain the phenotype...
I propose to do a new fluidgm on the time course 1(same samples as transcriptom samples) and not on time course 2 as the previous analysis to avoid difference and see more a effect on time course diversity

SNP-promoter diversity : make a short coherent list of genes to link it with domestication process?**



## Functional characterization of three AP2 genes

Three AP2 genes were selected because they are differentially expressed across species and their mutant were publicly available.

Those genes are *ERF48*, *ERF142* and *PLT8*.

![Expression of selected AP2 genes](../fig/fig-ap2-selected-expr.svg)
![Phenotype of Mutants selected AP2 genes](../fig/mutantAP2-boxplots.pdf)

But, it seems difficult ti finish with this point as :
ERF48 : present in cluster 2 and expression highly correlated with phenotype...but no phenotype observed in the Ox mutant
ERF142 : clear effect on panicle phenotype but in any modelisationPLT8 : effect on pbNb and SpNB..present in Cluster4..


# Discussion

Inflorescence branching in Solanaceae is modulated by subtle modifications of transcriptional programs [@Lemmon_Evolutioninflorescencediversity_2016].

# Materials and Methods

## Plant Material and Growth Conditions

For Phenotyping, 20 accessions ***(Supplementary Dataset 1)*** of *O. sativa*, *O. glaberrima*, *O. barthii* , plants were grown in the field of CIAT in Cali (Colombia) and in the greenhouse of IRD at Montpellier (France) for *O. rufipogon*.
At panicle maturity, the 3 main panicle from 3 plants per accession per repeat were collected (i.e 18 panicles/ accession). Each panicle was spread out and fixed on a white paper by tape. Panicle structure and seed number were analyses using P-TRAP software (Al-Tam et al., 2013). The quantified panicle traits included rachis length (RL in cm), total panicle length (PanL in cm), primary branch number (PbN), primary branch average length (PbL in cm), primary branch internode average length (PbintL in cm), secondary branch number (SbN), secondary branch average length (SBL in cm), secondary branch internode average length (SbintL in cm) and spikelet number (SpN). (Need a figure??). Phenotype description and statistical analysis of variance were performed by using functions in R software.
For expression analysis, plants were growth in greenhouse in Montpellier ...7 weeks in LD days and the shift to SD. 9 panicle mature per accessions have been phenotyped in Montpellier : table (phenotype_paniculesequeneced corrected.xls on folder Accession phenotype)

## Tissue collection and RNA sequencing
For expression analysis, around 15 panicles from each accession (IR64 for O. sativa indica, Niponbarre for O. sativa japonica, W1654 for O. rufipogon, Tog5681 for O. Glaberima, and B88 for O. Barthii) were collected from 4 days to 15 days after induction in 3 biological replicates for each repetition (kinetic 1 and 2).
Panicles were collected at 4 different morphological stages: stage 1, inflorescence (rachis ) meristem stage;  stage 2, branching stage (i.e, initiation of primary branches, panicles with elongated primary and higher order branch development); stage 3, spikelet differentiation; stage4, Floret differentiation and early flower organ development.


Total RNAs (including small RNAs) were extracted using RNeasy Plant Mini Kit with RLT and RWT buffers (Qiaegn, France). DNAse treatments were performed using the RNAeasy-free DNase set (Qiagen, France).

s

## Gene IDs and Names

Gene names were retrieved from [FunRiceGenes](https://funricegenes.github.io/) [@doi:10.1093/gigascience/gix119]

# Supplementary Material

- Supplementary Table 1: Panicle traits for 91 rice accessions.
- Supplementary Table 2: Independent measurements of panicle traits for 5 selected rice accessions.

# References

<div id="refs"></div>
