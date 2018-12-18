
## Materials and Methods

### Plant material and growth conditions

Panicle traits were measured in 91 accessions of *O. rufipgon*, *O. sativa*, *O. glaberrima*, *O. barthii* (Table S1 **Plantinfo**; **needs tidy**). Plants were grown in Cali, Colombia and Montpellier, France **move grow location to table S1**.
At panicle maturity, we collected the three main panicles from three plants per accession, per replicate (i.e. 18 panicles/ accession).
We used five accessions for expression analysis: *O. sativa japonica*  Nipponbare, *O sativa indica* IR64, *O. rufipogon* W1654, *O. glaberrima* Tog5681 and *O. barthii* B88.
These accessions were grown in a greenhouse in Montpellier, France, in June 2014, under long day conditions (14h light:10h dark).
After 6 to 8 weeks they were transferred to short day conditions (10h light:10h dark **is this right? 10 + 10 =/= 24**) to induce flowering.
To confirm panicle phenotypes in the growth conditions used for RNAseq, we evaluated panicle traits for 9 panicles from each accession that was grown in the greenhouse, under the same growth conditions.
The *crl5* and *smos1-3* mutants [@ayaNovelAP2TypeTranscription2014; @kitomiAuxinResponsiveAP22011] were grown in a greenhouse in Montpellier, France, in October 2017 under short day conditions (11h light:10h dark **is this right? 11 + 10 =/= 24**).
At least 18 panicles were used for panicle phenotyping.
All greenhouse plants were grown at 30°C with 80% relative humidity.
For phenotyping analyses, each panicle was spread out and fixed on a white paper using adhesive tape.
Panicle structure and seed number were analyses using `P-TRAP` software [@al-tamPTRAPPanicleTrait2013].

### Tissue collection and RNA sequencing

For expression analysis, we collected three biological replicates of approximately 15 immature panicles each from at least 10 plants per accession, per stage, from 4 days to 15 days after floral induction.
For sample collection, leaves surrounding the young panicle were removed by hand and the reproductive tissue was cut with a sharp blade under a (**brand, model**) stereomicroscope.
The reproductive tissues were frozen in liquid nitrogen, and total RNA including small RNA were extracted using RNeasy Plant Mini Kit with RLT and RWT buffers (QIAGEN, Germany).
DNase treatments were performed using the RNAeasy-free DNase set (QIAGEN, Germany).
RNA integrity numbers of the extracted RNA, measured using a 2100 Bioanalyzer (Agilent, U.S.A.), varied from 8.6–10.
Stage and meristem specificity were validated with quantitative real-time RT-PCR (qPCR) using stage-specific marker genes.
Primer sequences are listed in Table S2 **supp-table-PrimerList** **needs tidy**.
For each replicate, 400 ng of total RNA from the indeterminate (IM) and determinate (DM) meristem stages was used for RNAseq library preparation with the TruSeq Stranded Total RNA with Ribo-Zero Plant kit (Illumina, U.S.A.).
After quantification with a 2100 Bioanalyzer, 125-base paired-end reads were generated on a HiSeq 2500 (Illumina, U.S.A.) at the GeT Platform (Toulouse, Montpellier).

### Data analysis

#### Reproducibility and data availability

We used reproducible practices for all data analysis.
The code we used to analyse the RNAseq data and panicle phenotype data is hosted at https://github.com/tomharrop/5acc.
We used `snakemake` [@kosterSnakemakeScalableBioinformatics2012] to arrange analysis steps into workflows and monitor dependencies, and `Singularity` [@kurtzerSingularityScientificContainers2017] to capture the computing environment.
The results can be exactly reproduced from the raw data with a single command using `snakemake` and `Singularity`.

> Otho, can you please make sure the qPCR analysis complies with these reproducibility standards?

Raw sequence data are hosted at the NCBI SRA ... (**tom to do**)

#### Read alignment and Gene expression quantification

### Gene set enrichment
GSEA: Subramanian 2005, https://doi.org/10.1073/pnas.0506580102

FGSEA: Sergushichev 2016, https://doi.org/10.1101/060012

The pvalue of GSEA are estimated empirically by permutation.

### RT-qPCR
First-strand cDNA was synthesized from 1μg of DNase-treated total RNA using the SuperScriptIII cDNA First-stand synthesis system (Invitrogen) as recommended by the manufacturer and was stored at -20°C. Large-scale RT-qPCR was used to study the expression by taking advantage of the Microfluidic Dynamic Array and according to the manufacturer’s instructions  (Fluidigm, San Diego, CA, USA).
Before performing real-time PCR, the sample mixture and assay mixture were prepared individually. 96 × 96 Dynamic Array Integrated Fluidic Circuit was loaded with cDNAs (sample inlets) and primer combinations (assay inlets) after specific target amplification (15 cycles) (STA) and exonuclease I treatment. A fast cycling protocol and EvaGreen (Bio-Rad) as dye was used on a BioMark machine. The experiment was performed at the CIRAD platform (Montpellier, France) following the workflow provided by the manufacturer. Three biological replicates were performed for each sample. Data were normalized using 4 genes (LOC_Os06g11170, LOC_Os06g48970, LOC_Os01g16970, LOC_Os03g61680).
Relative gene expression were determined using the….methods

> TH: Otho, are you referring to the $2^{- \Delta \Delta C_{T}}$ method [@livakAnalysisRelativeGene2001]?

