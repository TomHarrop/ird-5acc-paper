
## Materials and Methods

### Plant material and growth conditions

Panicle morphological traits were measured in 91 accessions of *O. rufipogon*, *O. sativa*, *O. glaberrima* and *O. barthii*, grown in Cali, Colombia and Montpellier, France (Supporting information Table S1).
At panicle maturity, we collected the three main panicles from three plants per accession, per replicate (*i.e.* 18 panicles per accession).
We used four accessions for expression analysis: *O. sativa* ssp. *indica* IR64, *O. rufipogon* W1654, *O. glaberrima* Tog5681 and *O. barthii* B88.
These accessions were grown in a greenhouse in Montpellier, France, in June 2014, under long day conditions (14h light:10h dark).
After 6 to 8 weeks they were transferred to short day conditions (11h light:13h dark) to induce flowering.
To confirm panicle phenotypes in the growth conditions used for RNAseq, we evaluated panicle traits for 9 panicles from each accession that was grown in the greenhouse, under the same growth conditions.
The *crl5* and *smos1-3* mutants [@ayaNovelAP2TypeTranscription2014; @kitomiAuxinResponsiveAP22011] were grown in a greenhouse in Montpellier, France, in October 2017 under short day conditions (11h light:13h dark).
At least 18 panicles were used for panicle phenotyping.
All greenhouse plants were grown at 28°C with 80% relative humidity.
For phenotyping analyses, each panicle was spread out and fixed on white paper using adhesive tape.
Panicles were photographed and the images were used for panicle structure and seed number analysis with `P-TRAP` software [@al-tamPTRAPPanicleTrait2013].

### Tissue collection and RNA sequencing

For expression analysis, we collected 15 immature panicles each from at least 10 plants per accession, per stage, collected from 4 days to 15 days after floral induction (*i.e.* transfer to short day conditions).
For sample collection, leaves surrounding the young panicle were removed by hand and the reproductive tissue was cut with a sharp blade under a Stemi 508 (Zeiss, Germany) stereo microscope to identify developmental stage.
The reproductive tissues were immediately frozen in liquid nitrogen, and total RNA including small RNA was extracted using the RNeasy Plant Mini kit with RLT and RWT buffers (QIAGEN, Germany).
DNase treatments were performed using the RNase-free DNase set (QIAGEN, Germany).
RNA integrity numbers of the extracted RNA, measured using a 2100 Bioanalyzer (Agilent, U.S.A.), were between 8.6 and 10.
Stage specificity was validated with quantitative real-time RT-PCR (qPCR) using stage-specific marker genes (Supporting information Table S2).
400 ng of total RNA was used for each sample for RNAseq library preparation with the TruSeq Stranded Total RNA with Ribo-Zero Plant kit (Illumina, U.S.A.).
After quantification with a 2100 Bioanalyzer, 125-base paired-end reads were generated on a HiSeq 2500 (Illumina, U.S.A.) by the GeT Platform (Toulouse, France).

### qPCR

cDNA was synthesized from 1 μg of DNase-treated total RNA using the SuperScript III First-Strand Synthesis System (Invitrogen, U.S.A.). 
A Biomark HD Microfluidic Dynamic Array (Fluidigm, U.S.A.) was used for large-scale qPCR.
A 96 × 96 Dynamic Array Integrated Fluidic Circuit (Fluidigm, U.S.A.) was loaded with cDNA and primer combinations after 15 cycles of specific target amplification and exonuclease I treatment.
A fast cycling protocol with EvaGreen dye (Bio-Rad Laboratories, U.S.A.) was used for amplification.
Three biological replicates were performed for each sample.
Data were normalized using 4 genes (*LOC_Os06g11170*, *LOC_Os06g48970*, *LOC_Os01g16970*, *LOC_Os03g61680*).
Gene expression relative to the normalisation factors was estimated using the $2^{- \Delta \Delta C_{T}}$ method without a calibrator sample [@livakAnalysisRelativeGene2001].
Primer sequences are listed in Supporting information Table S2.

### Data analysis

We trimmed reads and removed adaptors with `cutadapt` [@martinCutadaptRemovesAdapter2011], before mapping to the MSU v7 annotation of the *Oryza sativa* ssp. *japonica* cv. Nipponbare reference genome [@ouyangTIGRRiceGenome2007] using `STAR` in 2-pass mode [@dobinSTARUltrafastUniversal2012].
To generate per-library gene expression cutoffs, we used the 95th percentile of reads that mapped to intergenic regions of the genome, as described previously [@harropGeneExpressionProfiling2016].
We used `DESeq2` [@loveModeratedEstimationFold2014] for differential expression analysis of genes that passed the cutoff.
We used annotations from the PlnTFDB v3.0 [@perez-rodriguezPlnTFDBUpdatedContent2010] and PlantTFDB v4.0 [@jinPlantTFDBCentralHub2017] to analyse expression of transcription factors.
Soft clustering of transcription factor genes was performed with `Mfuzz` [@kumarMfuzzSoftwarePackage2007], and enrichment of transcription factor family genes was tested with the GSEA method using the `FGSEA` package [@sergushichevAlgorithmFastPreranked2016; @subramanianGeneSetEnrichment2005]. 

### Reproducibility and data availability

Raw sequence data are hosted at the National Center for Biotechnology Information Sequence Read Archive (NCBI SRA) under accession [PRJNA518559](http://www.ncbi.nlm.nih.gov/bioproject/518559).
The code we used to analyse the RNAseq data and panicle phenotype data is hosted at https://github.com/tomharrop/5acc, and the code for qPCR analysis is at https://github.com/othomantegazza/5acc-qpcr.We used `snakemake` [@kosterSnakemakeScalableBioinformatics2012] to arrange analysis steps into workflows and monitor dependencies, and `Singularity` [@kurtzerSingularityScientificContainers2017] to capture the computing environment.
The final results and all intermediate steps can be exactly reproduced from the raw data with a single command using `snakemake` and `Singularity`.
The source for this manuscript is hosted at https://github.com/tomharrop/ird-5acc-paper.
