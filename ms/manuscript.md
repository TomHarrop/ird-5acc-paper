---
title: "Temporary - Determination of Inflorescence Architecture in Wild and Domesticated Rices"
bibliography: ../bib/references.bib
---

## Results

### Panicule traits that influence spikelet number

Besides the reference *Oryza sativa nipponbare*, we have selected 4 other rice accessions, one per species: B88 [*Oryza barthii*], Tog5681 [*Oryza glaberrima*], W1654 [*Oryza Rufipogon*], IR64 [*Oryza sativa indica*], because they are easily grown in greenhouse, because they have sequenced or nearly sequenced genomes and most important because they display representative phenotypic traits.

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

We asked: what are the panicule traits that most influence spikelet number in all species? Thus we considered spikelet number as an outcome and all the other phenotypic traits as potential predictor of a linear model and we have rephrased our question in stricter statistical terms as what is the minimum optimal number of predictor needed to explain and predict the outcome?

We have answered that question with linear model subset selection [@miller2002subset, @lumey2017leapspackage]. This algorithm searches for the combination of predictor that it is likely to provide the most accurate prediction of the outcome. The suggested combination is: primary branch number, secondary branch number and primary branch length, with the effect of the last term being small.

Thus, we have fitted a linear model on primary and secondary branch number. The coefficients of that linear model successfully predict spikelet number in an independent dataset. This independent dataset reports the same parameters for 100 rice accession from the 4 species taken in consideration (*barthii*, *glaberrima*, *rufipogon* and *sativa*). Since in 100 accessions we are able to predict spikelet number from only three parameters we show that the rules that determine the spikelet number within panicules of different rice species have the same basis (*fig-spn-prediction*).

![Spikelet number is predicted by primary and secondary branch number](../fig/fig-spn-prediction.svg)

The two first principal components of this dataset explain 28.6% and 16.3% of variance

### PCA: the 5^th^ component splits developmental stages

We have used Principal Component Analysis (PCA) in order to explore and detect patterns of gene expression across different species and stages. The first two component split clades (Asian vs African Rice in the 1^st^ and 2^nd^ components); the next two components split species within the same clade (Rufipogon vs. Sativa in the 3^rd^; Barthii vs. Glaberrima in the 4^th^). Indeed species are the first source of variation, but this source of variation could be confounded with mapping bias.

Further components can be explained by differences between the two developmental stages: PBM and SM.  Above all, the 5^th^ component splits PBM and SM in all species (*fig-pc5*).

![fig-pc5: The 5^th^ component splits developmental stages ](../fig/fig-pc5.svg)

### PCA: master regulators of development are enriched in the loadings of the 5^th^ component

The 5^th^ principal component must be driven by gene expression patterns that identify the **primary branch meristem** and the **spikelet meristem** in all or most species. Are those pattern biologically significant? We have described those pattern by using the loading of each gene on the 5^th^ component as rankings and checking if they enrich mapman functional categories [@thimm2004mapman].

The rankings of the 5^th^ component enrich multiple mapman bins (*table-pc5-mapman*). Many of the enriched mapman bins are related to inflorescence development.

![table-pc5-mapman](../tables/table-pc5-mapman.svg)

### Functional characterization of three AP2 genes

Three AP2 genes were selected because they are differentially expressed across species and their mutant were publicly available. *EXTEND!*

Those genes are *ERF48*, *ERF142* and *PLT8* (*fig-ap2-selected-expr*).

![Expression of selected AP2 genes](../fig/fig-ap2-selected-expr.svg)

## Discussion

Inflorescence branching in Solanaceae is modulated by subtle modifications of transcriptional programs [@Lemmon_Evolutioninflorescencediversity_2016].

## References

<div id="refs"></div>
