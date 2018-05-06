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

## Diversity of panicle architecture in the Asian and African Rices

Phenotyping of the 20 accessions in the 4/5 species : PCA and correlation analysis
[fig-PCA diversity: Panicle traits architecturehcitectire](../fig/Test FigPCA.pdf)
2 conclusions :
 - Panicle traits that influences panicle architecture Diversity and spikelet number
 - Panicle traits correlated with Sp nb
 - Obersvation of panicle architecture similarity between the 2 domestications process

> **Presentation of the transcriptomics parts - Phenotype of the 5 accessions  : Let's see how to introduce this/ it will depend also if we include or not the analysis by graph from the modelisation group**


Besides the reference *Oryza sativa nipponbare*, we have selected 5 other rice accessions, one per species: B88 [*Oryza barthii*], Tog5681 [*Oryza glaberrima*], W1654 [*Oryza Rufipogon*], Niponbarre [*Oryza sativa japonica], IR64 [*Oryza sativa indica*], because they are easily grown in greenhouse, because they have sequenced or nearly sequenced genomes and most important because they display representative phenotypic traits.

On those accessions, we have measured systematically 9 panicle traits [@faroq2013p]. These traits are:

- Rachis Length,
- Primary Branch Numbers,
- Primary Branch Average Length,
- Primary Branch Internode Average length,
- Secondary Branch Numbers,
- Secondary Branch Numbers,
- Secondary Branch Average Length,
- Tertiary Branch Numbers,
- Spikelet Numbers.

Panicule architecture influences spikelet number (***This is kind of obvious, but do we have a good reference?***),  which is strongly increased during domestication. Different species might use different strategies to lodge more spikelets on one panicule. For example sativa indica relies more on secondary branches, while sativa japonica relies on both on primary and secondary branches and glaberrima relies more on primary branches (*fig-branches-boxplot*).

![fig-branches-boxplot: Different species use different panicule architecture to lodge more spikelets](../fig/fig-branches-boxplot.svg)
also the previous [fig-phenotype: panicle traits in the 5 accessions]fig (../fig/traits_corrected.pdf)

We asked: what are the panicule traits that most influence spikelet number in all species? Thus we considered spikelet number as an outcome and all the other phenotypic traits as potential predictor of a linear model and we have rephrased our question in stricter statistical terms as what is the minimum optimal number of predictor needed to explain and predict the outcome?

We have answered that question with linear model subset selection [@miller2002subset, @lumey2017leapspackage]. This algorithm searches for the combination of predictor that it is likely to provide the most accurate prediction of the outcome. The suggested combination is: primary branch number, secondary branch number and primary branch length, with the effect of the last term being small.

Thus, we have fitted a linear model on primary and secondary branch number. The coefficients of that linear model successfully predict spikelet number in an independent dataset. This independent dataset reports the same parameters for 100 rice accession from the 4 species taken in consideration (*barthii*, *glaberrima*, *rufipogon* and *sativa*). Since in 100 accessions we are able to predict spikelet number from only three parameters we show that the rules that determine the spikelet number within panicules of different rice species have the same basis (*fig-spn-prediction*).

![Spikelet number is predicted by primary and secondary branch number](../fig/fig-spn-prediction.svg)

The two first principal components of this dataset explain 28.6% and 16.3% of variance


Presentation of the stage/time course, checking some genes by Q-PCR (perhaps other genes necessary which could be include in the next fluidgm)

> **What are the genes/patwhay involved in the meristem state in rice inflorescence Undeterminate vs determinate**

## PCA: the 5^th^ component splits developmental stages

We have used Principal Component Analysis (PCA) in order to explore and detect patterns of gene expression across different species and stages. The first two component split clades (Asian vs African Rice in the 1^st^ and 2^nd^ components); the next two components split species within the same clade (Rufipogon vs. Sativa in the 3^rd^; Barthii vs. Glaberrima in the 4^th^). This indicates that species are the first source of variation of genes expression, but this source of variation could be confounded with mapping bias.

Further components can be explained by differences between the two developmental stages: PBM and SM.  Above all, the 5^th^ component splits PBM and SM in all species (*fig-pc5*).

![fig-pc5: The 5^th^ component splits developmental stages ](../fig/fig-pc5.svg)

## The 5^th^ principal component is driven by master regulators of development

The 5^th^ principal component must be driven by gene expression patterns that identify the **primary branch meristem/axillary meristem** and the **spikelet meristem** in all or most species. Are those pattern biologically significant? We have described those pattern by using the loading of each gene on the 5^th^ component as rankings and checking if they enrich mapman functional categories [@thimm2004mapman].

The rankings of the 5^th^ component enrich multiple mapman bins (*table-pc5-mapman*). Many of the enriched mapman bins are related to inflorescence development.

![table-pc5-mapman](../tables/table-pc5-mapman.svg)

### Within the 5^th^ component, *indica* IR64 is an outlier

![](../fig/fig-genes-of-pc5.svg)

The other species cluster together randomly

## The genes that drive PC5 might determine spikelet development and its divergence

Are the differences determined in the PBM (top pc5 genes) or in the SM (bottom pc5 genes)?

### PBM - Top PC5 - Indica diff

![](../fig/fig-TMP-pc5-top200-cl1.jpeg)

OsLHY, RDD1, **Circadian Clock** behaves different in indica, but also in other asian species.

![](../fig/fig-TMP-pc5-top200-cl4.jpeg)

OsRBCS4, **Photosynthesis** behaves differently in indica, different development timing? Also: OsPORA,

OsCML16,  **enhances root growth** flat in indica and japonica

OsABA8ox1 **ABA hydrolase** regulates growth, behaves opposite in indica

OsDREB1E, AP2-like **Osmotic stress / Downstream of Gigantea?**,
OSISAP1, **resistance to dehyd.**,

OsNAC9, **Root structure**,

OsACO1, **Internode elongation**

OsRAV2, **regulated by salt?**

OsBIERF4, **Something stress again??**

AP37, Stress again, **programmed cell death**

OsMYL1, bHLH something **cell identity**.

Other cluster, still ***Indica*** peculiar:

![](../fig/fig-TMP-pc5-top200-cl5.jpeg)

OsDOS: CCCH **Senescence** Something AP2s network???

IAA6|OsIAA6: **auxin throuh PIN1?** tiller outgrow + **stress**

OsPIP1: **stress** auquaporin

GLN1;2|GS1;2|OsGS1;2: **high temperature** tolerance

OsGASR1: **Differentiation of panicle**

cZOGT1: **Citokinin related**

OsMADS16: **mads B-class**

OsSAP11: **citokinin - cell elongation**

### PBM - Top PC5 - Indica same

![](../fig/fig-TMP-pc5-top200-cl3.jpeg)

Interesting genes highlighed : AP2/ALOG on one side and MADS box genes in the other side of the FigPCA. There is also genes related to hormone.

Interestingly, some genes seems to have a a variant pattern in indica (some AP2 and MADS box) : these genes are present on cluster3 and 4
Maybe make a bridge with phenotype (high and low branching) that could be explained by these genes.

### Cluster analysis

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

Three AP2 genes were selected because they are differentially expressed across species and their mutant were publicly available. *EXTEND!*

Those genes are *ERF48*, *ERF142* and *PLT8* (*fig-ap2-selected-expr*).

![Expression of selected AP2 genes](../fig/fig-ap2-selected-expr.svg)
![Phenotype of Mutants selected AP2 genes](../fig/mutantAP2-boxplots.pdf)

But, it seems difficult ti finish with this point as :
ERF48 : present in cluster 2 and expression highly correlated with phenotype...but no phenotype observed in the Ox mutant
ERF142 : clear effect on panicule phenotype but in any modelisationPLT8 : effect on pbNb and SpNB..present in Cluster4..


# Discussion

Inflorescence branching in Solanaceae is modulated by subtle modifications of transcriptional programs [@Lemmon_Evolutioninflorescencediversity_2016].

# References

<div id="refs"></div>
