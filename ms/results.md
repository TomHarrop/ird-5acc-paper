## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice, followed by whole-transcriptome sequencing (RNASeq) [and analysis of single-gene mutants with defects in panicle branching].
**expand slightly**

### Parallel changes in panicle architecture in domesticated accessions

To measure the diversity of panicle architecture, we phenotyped 93 accessions of wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*), using `P-TRAP` for automated measurement of traits [@al-tamPTRAPPanicleTrait2013] (Figure 1A-D - **phenotype-pca**).
There was a stronger correlation between secondary branch number and spikelet number in the domesticated accessions than in the wild accessions (Figure S1 **correlation-pbn-spn**).
The first principal component (PC1) in the phenotyping data accounts for 46.5% of variability, and separates domesticated and wild accessions but not Asian and African accessions (Figure 1E **phenotype-pca**).
The sum of contributions of the top four components accounted for more than 87% of the total variance, but components other than PC1 do not separate panicles from different accessions (Figure S2 **suppl-phenotype-pca-all-pc**), and PC1 is the only component that splits the accessions by domestication status.
Spikelet number, secondary branch number and primary branch number have the highest loadings on PC1, suggesting that these factors, which are all related to axillary meristem formation and fate transition, are the main differences between panicles from wild and domesticated accessions.


### Differences in gene expression profiles between wild and domesticated panicles

To investigate gene expression associated with variation in panicle architecture, we chose a single accession for each of *O. rufipogon* (W1654), *O. sativa japonica* (Nipponbare), *O. sativa indica* (IR64), *O. barthii* (B88) and *O. glaberrima* (Tog5681) for RNA sequencing (RNAseq).
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S3 **suppl-phenotype-all-varieties**), except for *O. sativa japonica* cv. Nipponbare, which we included as the reference accession for *O. sativa japonica*.
> OM: Helene, maybe you can extend a bit this part on selecting accessions
We performed further panicle phenotyping on 9 panicles from 3 plants from the accessions chosen for RNAseq.
The three domesticated accessions produce more spikelets than their wild relatives (Figure S4 **phenotyping-mpl**).
**Model of panicles goes here.**

To measure whole-transcriptome gene expression, we collected immature panicles from each accession at four developmental stages: rachis meristem (RM); branch meristem (BM), composed of immature panicle displaying primary branch initiation, elongation of primary branches and axillary meristem initiation; spikelet meristem (SM), which includes panicles with spikelet meristem, but not floret meristem; and flower meristem (FM) (**S5a qpcr-confirms-sampling**). 
We extracted RNA from pooled meristems at each stage, and measured expression of *LHS1*, *G1L5* (*TAWAWA*), *FZP*, *LAX1* and *MADS14* with quantitative RT-PCR (qRT-PCR) to confirm staging (Figure S5b **qpcr-confirms-sampling**).
Because branching complexity is related to branch meristem establishment and meristem fate transition, we chose the BM and SM stages for RNAseq.
cDNA libraries for sequencing were constucted from rRNA-depleted RNA from these stages for all five accessions.
We mapped reads against the MSU Release 7.0 of the annotation of the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007].
We obtained an average of more than 20 million uniquely mapped reads within exons for each accession, including African rice species (Table S1 **mapping-statistics**).
> **is it possible to add more information about the % of rRNA is the libaraies..as these libararies are not the common used?**
> No, it's not interesting for this paper. rRNA depletion is the standard alternative to polyA enrichment.
Pairwise distances between libraries were lower for samples from the same stage, accession and continent, in that order (Figure S6 **distance-heatmap**), indicating that transcriptomic changes during domestication are subtle compared to differences between species.

We used differential expression (DE) analysis to identify genes that were up- or down-regulated between stages across all accessions.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1 (Table S2 **DE-genes-stages**).
This list included **TFs x, y and z** [e.g. pick from *FZP*, *LHS1*, *LAX1*, *PAP2*, *MFO1*, ...], which control inflorescence architecture in rice [**add citations for chosen genes**], **indicating that our approach identifies genes that control branching.**

## AP2/EREBP-like transcription factors are differentially expressed between stages and associated with domesticated accessions

We used PCA on transformed raw counts to investigate the patterns of variation in the transcriptomes. 
The first four PCs split different combinations of rice species (Figure 2 **transcriptome-pca**). **expand?**
> I don't think so, these PCs are not interesting
These components may relate to species-specific differences unrelated to panicle architecture, or mapping biases caused by mapping all libraries against the *O. sativa japonica* reference.
In contrast, PC5 separates developmental stages across all five accessions except *O. sativa indica* (Figure 2 **transcriptome-pca**). **expand?**
> What else is there to say about PC5?
The family of APETALA2 and ethylene-responsive element binding protein (AP2/EREBP)-like genes are enriched at the extremes of PC5, suggesting that they contribute to differences between BM and SM (Figure 3 **HB-AP2-heatmap**).
Homeobox, MADS, NAC and SBP genes are also enriched in PC5 (Figure 3 **HB-AP2-heatmap**; Figure S7 **NAC-MADS-SPL-heatmap**).

Because of the prominence of transcription factor (TF) genes and AP2/EREBP genes in PC5 and in DE genes between stages (Figure 3 **HB-AP2-heatmap**; Table S2 **DE-genes-stages**), we used soft clustering of log~2~-fold change values (L~2~FCs) between BM and SM to find common patterns of expression of the subset of annotated TF genes that were expressed in our RNAseq dataset.
There were six common patterns, three of which had the highest L2FC, indicating higher expression in SM than in BM, in *O. sativa indica*.
These three clusters were all correlated with secondary branch number and spikelet number (Figure 4 **cluster-phenotype-corr**).
Cluster 5, which had the highest core L~2~FC in *O. sativa indica* and the highest correlations with secondary branch number and spikelet number, had an enrichment of AP2/EREBP genes (two-tailed hypergeometric test; *p*~adj~ == **x**).
One possibility is that genes in this cluster specify SM, and their delayed repression in *O. sativa indica* results in more spikelets and higher-order branches in this species.
Most of the genes in cluster 5 have negative L~2~FCs in the other accessions, indicating repression of these genes in the SM stage, and L~2~FCs close to zero in *O. sativa indica*, consistent with lack of repression in that accession (Figure S8 **cluster-5-details**).
**Mention other clusters**.

To find TF genes associated with changes in panicle architecture during domestication, we tested the stage × accession interaction for African and Asian accessions separately at an FDR of 0.1 (Table S3 **DE-genes-interaction**).
For Asian accessions (*O. rufipogon* and *O. sativa indica*), there was a significant interaction for 85 genes, including 12 AP2/EREBP-like genes.
In African accessions, the stage × accession interaction was significant for 50 genes, including 8 AP2/EREBP-like genes (**check ap2 numbers**).
*INDETERMINATE SPIKELET 1* (*IDS1*), which controls inflorescence architecture (**ref**), was DE in both comparisons, consistent with its importance in rice domestication (**check if IDS1 paper mentions domestication**).
The other genes in these lists are candidate targets of artificial selection for changes in panicle architecture.
Expression of the 10 genes that appear in both comparisons, including *IDS1* and the AP2/EREBP-like gene *ERF74*, may have evolved in parallel in the separate domestication of Asian and African rice.

## AP2/EREBP mutants have defects in panicle branching

Next, we wanted to test whether loss-of-function mutants had phenotypes consistent with our hypothesis that changes in expression of AP2/EREBP-like genes control panicle architecture.
Mutant generation in rice takes up to **x** months/years (**reference**), so we chose to start by phenotyping panicles in mutants that were publicly available.
A **knockout?** mutant of *PLETHORA 8* (*PLT8*), an AP2-like gene characterised in (**ref**), produces a shorter rachis with fewer primary branches than the background accession (Figure 5 **panicle-mutants**), consistent with its reported peak in expression in rachis meristem tissue in *O. sativa japonica* cv. Nipponbare [Figure S9 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
*ERF142* and *DLT* mutants, reported in (**ref**), produce fewer primary and secondary branches and fewer spikelets (Figure 5).
>**in discussion**: We detected all three genes at both stages of all five accessions in our RNAseq dataset, but none of them were differentially expressed. Although this suggests that they were not targets of domestication, the phenotypes support a role for AP2/EREBP-like genes in panicle architecture.
> Could we say these genes are involved in branching but not necessarily domestication - that's why they didn't show up in our analysis?
