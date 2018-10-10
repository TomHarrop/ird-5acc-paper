## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice, followed by whole-transcriptome sequencing (RNASeq) [and analysis of single-gene mutants with defects in panicle branching].
**expand slightly**

### Convergence of panicle architecture between the 2 domestications process. 

The panicle architecture was phenotyped in 93 accessions of wild Asian rice (*Oryza rufipogon*), domesticated Asian rice (*Oryza sativa*), wild African rice (*Oryza barthii*) and domesticated African rice (*Oryza glaberrima*)(**Supp-Table01-Plantinfo**). Based on panicle images and `P-TRAP` software for automated measurement of panicle traits [@al-tamPTRAPPanicleTrait2013] (Figure 1A-D - **phenotype-pca**), a total of 10 panicle traits were scored (**Supp-Table02-PanicleTraitsPhenotyping**). Principal components analysis (PCA) was used to elucidate the variation of panicle architecture between rice accesssions. These accessions could be grouped depending of several characteristics: species, geographical origin (Asia/Africa) and type (wild/domesticated) . 
The first principal component (PC1) in the phenotyping data accounts for 46.5% of variability The sum of contributions of the top four components accounted for more than 87% of the total variance. The contributions of PC2, PC3 and PC4 are 18,3%, 12,3% and 10,3% respectively. The most loadings of the PC1 are Spikelet number (SpN), secondary branch number (SBN) and primary branch number (PBN). For PC2, the main loadings are Primary branch Number (PBN), and the lenght of Primaryt branch (PBL) and Secondary Branch (SBL). The PC3 is supported by Rachis lengtn (RL), Tertiary Branch number (TBN) and the internode length along the Secondary branch (SBIL) (**SF1**). 
Projection of the accessions in all four component was used to analyse the relationships between panicle architecture and the differents characteritics of these accessions (species/geographical orifin and type). The first components is the only components able to split the accessions grouped by the type characteristics (wild/domesticated) (Figure 1E **phenotype-pca**). The panicle archtiecture is not able to distinguish clearly the accessions based on the species and/or geographical origin chatcateristics (**SF1** suppl-phenotype-pca-all-pc). 
>Be sure it is enough clear based on our figures.

Correlations among most loadings of PC1 (SpN, PBN, SBN) were analyzed separatly in the 4 species (Figure S1 **correlation-pbn-spn**) to evaluate. WE obersved a strong positive correlation between the Spn and SBN in the domesticated accessions that in wild accessions. In contrast, a higher positive correlation is observed between SpN and PBN in the wild accessions than in domesticated accessions. 

>please add the correlation value in each panel of the figure correlation-pbn-sbn.

CCL of this part : 
- PC1 (PBN, SBN, SpN) (PBN, SBN, SpN) (which are relative to axillary meristem) are the main contributors of panicle diversity
- only the relationships between Type and panicle architecture is able to groupe the accession. suggesting a similar branching phenotype beteen the wild and the domesticated accessions. 
- the differnet correlation between SpN, SBN and PBN between the wild and domesticatexd suggest different control of branching pattern.


### Panicle transcriptome Variation (Differences in gene expression profiles between wild and domesticated panicles : **remove this, we don't show it**)

To investigate gene expression associated with variation in panicle architecture, we chose a single accession for each of *O. rufipogon* (W1654), *O. sativa japonica* (Nipponbare), *O. sativa indica* (IR64), *O. barthii* (B88) and *O. glaberrima* (Tog5681) for RNA sequencing (RNAseq). 
Based on the extensive phenotyping described above, the chosen accessions are consistent with species-wide patterns of panicle architecture (Figure S3 **suppl-phenotype-all-varieties = suppl-fig-MAYBE-pca-on-pheno-BOXPLOT-LABEL in dropbox have to be finalize**), except for *O. sativa japonica* cv. Nipponbare, which we included as the reference accession for *O. sativa japonica*.
> OM: Helene, maybe you can extend a bit this part on selecting accessions  (** i'm still waiting the mapping % of rufipogon against IR64 and NIP...and i was thinking that could be interesting to have a SNP diversity analysis also...i'm waiting Francois sabot and Christine to see how to manage this, it will be in suppdata**)

To measure whole-transcriptome gene expression related to panicle branching diversity, we collected immature panicles from each accessions at different developmental stages : stage 1, rachis meristem (RM); stage2, branch meristem (BM) composed by immature panicle displayong primary branch initiation and elongation of these primary branch with axxilary meristem initiation; ; Stage3, spikelet meristem (SpM) which is a mix of early Spikelet meristem differentiation until floret meristem; and stage4, flower meristem (FM) stages (Figure S5a **qpcr-confirms-sampling** **Finalize this figure**).
We extracted RNA from pooled meristems of the same stage for the 5 accessions in three biological replicates. Staging were validated with quantiative RT-PCR using stage specific marker genes : *LHS1* ; *G1L5* FZP, LAX1, OsMADS14, (*TAWAWA*) (Figure S5b **qpcr-confirms-sampling**) (**add description**). Panicle traits were also investigated for 9 panicles from 3 plants of the accessions chosen for RNAseq and grow in the same conditions. The three domesticated accessions produce more spikelets than their wild relatives **add more details in the description**, with variation in Sp, Sb, pb nb between the 5 species** (Figure S4 **phenotyping-mpl** **Supp-table03-panicletrairsPhenotypingPlants Sequenced**).
**Model of panicles goes here : HA, i will prepared something and see if it is interesting to add it or not**. 

As the branching pattern is related to the stage of branch meristem establishment and meristem fate transition, we choose to used the RNA from the BM (Stage 2) and SpM (Stage3) stages for performing libraries using the TrueSeq Stranded Total RNA with Ribo-Zero Plant prep Kit and performed illumina sequencing. 
We mapped reads against the MSU Release 7.0 of the annotation of the *O. sativa japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007], and we obtained an average of more than 20 million uniquely mapped reads within exons for each accession per biological replicate, including African rice species (Table S1 **mapping-statistics**). **is it possible to add more information about the % of rRNA is the libaraies..as these libararies are not the common used?**
Pairwise distances between libraries were lower for samples from the same stage, accession and continent, in that order (Figure S6 **distance-heatmap**), indicating that transcriptomic changes during domestication are subtle compared to differences between species.

To identify patterns variation in transcriptomes of the samples, we used PCA analyses **explained more which raw data were used; raw counts + quantification or not....**
The first four PCs split different combinations of rice species (Figure 2 **transcriptome-pca**) **describe more variance....**
These components may relate to species-specific differences unrelated to panicle architecture, or mapping biases caused by mapping all libraries against the *O. sativa japonica* reference. **describe a  more the 4 components**
In contrast, PC5 separates developmental stages across all five accessions except *O. sativa indica* (Figure 2 **transcriptome-pca**). **moderate the less separation in indica**

Then, we used differential expression (DE) analysis to identify genes that were up- or down-regulated between stages across all accessions.
193 genes were at least 1.5-fold DE between stages in all species at a false-discovery rate of 0.1 (Table S2 **DE-genes-stages**).
This list included **TFs x, y and z** [e.g. pick from *FZP*, *LHS1*, *LAX1*, *PAP2*, *MFO1*, ...], which control inflorescence architecture in rice [**add citations for chosen genes**].
> Shall we discuss the fact that many more genes are differentially expressed in African species? (also supported by the PCA) Or would this make the paper overcomplicated? Shall we discuss differential expression in the single species at all?
> descroibe more the resulsts

CCL : **?**

## AP2/EREBP-like transcription factors are differentially expressed between stages and associated with domesticated accessions
>make the link with the previous part from all genes to TF? supp-fig-05-tfs-of-pc5

The family of APETALA2 and ethylene-responsive element binding protein (AP2/EREBP)-like genes are enriched at the extremes of PC5, suggesting that they contribute to differences between BM and SM (Figure 3 **HB-AP2-heatmap**).
Homeobox, MADS, NAC and SBP genes are also enriched in PC5 (Figure 3 **HB-AP2-heatmap**; Figure S7 **NAC-MADS-SPL-heatmap**).
> what is the link with domesticated here?
> do we refer to LMD data/AP2-HB in these part?

## Transcriptome Association with panicle complexity
Because of the prominence of transcription factor (TF) **we don't show really this prominence..or give a statistic value or comparison between all genne/TF** genes and AP2/EREBP genes in PC5 and in DE genes between stages (Figure 3 **HB-AP2-heatmap**; Table S2 **DE-genes-stages**), we used soft clustering of log~2~-fold change values (L~2~FCs) between BM and SM to find common patterns of expression of the subset of annotated TF genes that were expressed in our RNAseq dataset.
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

## Differences in AP2/EREBP promoters are associated with changes in expression pattern in domesticated species

> **removed this part** **In progress?**

## AP2/EREBP mutants have defects in panicle branching

Next, we wanted to test whether loss-of-function mutants had phenotypes consistent with our hypothesis that changes in expression of AP2/EREBP-like genes control panicle architecture.
Mutant generation in rice takes up to **x** months/years (**reference**), so we chose to start by phenotyping panicles in mutants that were publicly available.
A **knockout?** mutant of *PLETHORA 8* (*PLT8*), an AP2-like gene characterised in (**ref**), produces a shorter rachis with fewer primary branches than the background accession (Figure 5 **panicle-mutants**), consistent with its reported peak in expression in rachis meristem tissue in *O. sativa japonica* cv. Nipponbare [Figure S9 **lmd-paper-ap2**; @harropGeneExpressionProfiling2016].
*ERF142* and *DLT* mutants, reported in (**ref**), produce fewer primary and secondary branches and fewer spikelets (Figure 5).
>**in discussion**: We detected all three genes at both stages of all five accessions in our RNAseq dataset, but none of them were differentially expressed. Although this suggests that they were not targets of domestication, the phenotypes support a role for AP2/EREBP-like genes in panicle architecture.
> Could we say these genes are involved in branching but not necessarily domestication - that's why they didn't show up in our analysis?
