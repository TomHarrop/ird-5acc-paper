---
title: "Temporary - Determination of Inflorescence Architecture in Wild and Domesticated Rices"
bibliography: ../bib/references.bib
---

## Results

### Seed number is controlled by how many independent pathways?

### The minimum number of pathways that control panicle architecture

Besides the reference *Oryza sativa nipponbare*, we have selected 4 other rice accessions, one per species: B88 [*Oryza barthii*], Tog5681 [*Oryza glaberrima*], W1654 [*Oryza Rufipogon*], IR64 [*Oryza sativa indica*]. Because they are easily grown in greenhouse, they have sequenced or nearly sequenced genomes and most important they display representative phenotypic traits.

On those accessions, we have measured systematically 9 panicle traits [@faroq2013p]. These traits are:

- Rachis Lengths,
- Primary Branch Numbers,
- Primary Branch Average Lengths,
- Primary Branch Internode Average lengths,
- Secondary Branch Numbers,
- Secondary Branch Numbers,
- Secondary Branch Average Lengths,
- Tertiary Branch Numbers,
- Spikelet Numbers.

Arguably, spikelet number is one of the main outcomes of panicle architecture and one of the main goals of domestication. Thus we consider spikelet number as an outcome and all the other phenotypic traits as potential predictor.

In order to infer the minimum number of independent pathway that determine spikelet number during panicle development, we asked how many panicle traits are necessary to predict spikelet number?

To identify those traits we have performed linear model selection.

[.... *more details on model selection here* ....]

The two first principal components of this dataset explain 28.6% and 16.3% of variance

### PCA: the 5^th^ component splits developmental stages

We have used Principal Component Analysis (PCA) in order to explore and detect patterns of gene expression across different species and stages. The first two component split clades (Asian vs African Rice in the 1^st^ and 2^nd^ components); the next two components split species within the same clade (Rufipogon vs. Sativa in the 3^rd^; Barthii vs. Glaberrima in the 4^th^). Indeed species are the first source of variation, but this source of variation could be confounded with mapping bias.

Further components display patterns that can be explained by differences between the two developmental stages: PBM and SM.  Above all, the 5^th^ component splits PBM and SM in all species (*fig-pc5*).

![fig-pc5](../fig/fig-pc5.svg)

### PCA: master regulators of development are enriched in the loadings of the 5^th^ component

The 5^th^ principal component must be driven by gene expression patterns that identify the **primary branch meristem** and the **spikelet meristem** in all or most species. Are those pattern biologically significant? We have described those pattern using the loading of each gene on the 5^th^ component as ranking and checking if those rankings enrich functional categories.  

The rankings of the 5^th^ component enrich multiple mapman bins (*table-pc5-mapman*). Many of the enriched mapman bins are related to inflorescence development

![table-pc5-mapman](../tables/table-pc5-mapman.svg)

## Discussion

Inflorescence branching in Solanaceae is modulated by subtle modifications of transcriptional programs [@Lemmon_Evolutioninflorescencediversity_2016].

## References

<div id="refs"></div>
