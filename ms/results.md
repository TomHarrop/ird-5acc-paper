## Results

To reveal gene expression patterns associated with diversity of inflorescence architecture, we used detailed phenotyping of panicles from domesticated and wild accessions of Asian and African rice, followed by whole-transcriptome sequencing (RNASeq) [and analysis of single-gene mutants with defects in panicle branching]. **expand slightly**

### Panicles from domesticated accessions produce more primary branches, secondary branches and spikelets

To measure the diversity of panicle architecture, we phenotyped **n** accessions of wild Asian rice (*Oryza rufipogon*; figure 1A), domesticated Asian rice (*Oryza sativa*; figure 1B), wild African rice (*Oryza barthii*; figure 1C) and domesticated African rice (*Oryza glaberrima*; figure 1D), using `P-TRAP` for automated measurement of traits [@al-tamPTRAPPanicleTrait2013]. The first principal component in the phenotyping data accounts for 46.5% of variability, and separates domesticated and wild accessions but not Asian and African accessions (figure 1E). Spikelet number, secondary branch number and primary branch number have the highest loadings on PC1, suggesting that these factors are the main differences between panicles from wild and domesticated accessions. Lower-order components do not separate panicles from different accessions (figure S1). **Correlations**.

### Differences in gene expression profiles between wild and domesticated panicles

To investigate gene expression associated with variation in panicle architecture, we chose a single accession for each of *O. rufipogon*, *O. sativa* japonica (Nipponbare), *O. sativa* indica (IR84), *O. barthii* (**?**) and *O. glaberrima* (**TOG...**). We performed further phenotyping on these accessions to ensure that they were consistent with species-wide patterns of panicle architecture (Figure S3).

## Figure legends

**Figure 1**. Panicles from domesticated accessions produce more branches and spikelets than panicles from wild accessions. We used spread panicles from **A** *O. rufipogon*, **B** *O. sativa*, **C** *O. barthii* and **D** *O. glaberrima* to measure panicle phenotypes with P-TRAP [@al-tamPTRAPPanicleTrait2013]. The first principal component (PC1) in the panicle phenotye data accounts for 46.5% of variability and separates wild and domesticated accessions (**E**), and spikelet number (SpN), secondary branch number (SBN) and primary branch number (PBN) had the highest loadings on PC1 ( **F**). RL: Rachis length; PBL: Primary branch length; PBIL: Primary branch internode length; SBL: Secondary branch length; SBIL: Secondary branch internode length; TBN: Tertiary branch number; PL: Panicle length.

**Figure S1**. Principal component analysis (PCA) of panicle phenotyping data showing components 1–4. PC1 accounts for 46.5% of variability and separates panicles from domesticated and wild accessions. The lower ordinates do not separate panicles by species.

**Figure S2** Correlation between the main traits that determine panicle phenotype. A. Primary Branch Number and Spikelet Number correlates only slightly mainly in wild species. B. Secondary Branch Number and Spikelet number
highly correlates in cultivated species, less in wild species. C.
Primary and Secondary Branches poorly correlates, suggesting that they
are controlled by different genetic mechanisms



## Materials and Methods

###Plant Material and Growth Conditions
Panicle traits diversity were investigated in accessions (**Supp-Table S1**) of wild and domesticated, Asian and African rice.  O. sativa, O. glaberrima, O. barthii  and o. plants were grown in Cali (Colombia) and in Montpellier (France) for O. rufipogon. 
At panicle maturity, the 3 main panicle from 3 plants per accession per repeat were collected (i.e. 18 panicles/ accession). Each panicle was spread out and fixed on a white paper by tape. Panicle structure and seed number were analyses using P-TRAP software (Al-Tam et al., 2013). The quantified panicle traits included rachis length (RL in cm), total panicle length (PanL in cm), primary branch number (PbN), primary branch average length (PbL in cm), primary branch internode average length (PbintL in cm), secondary branch number (SbN), secondary branch average length (SBL in cm), secondary branch internode average length (SbintL in cm) and spikelet number (SpN). Phenotype description and statistical analysis of variance were performed by using functions in R software.
For expression and sequencing analysis, 5 accessions have been used (Niponbarre for O. sativa japonica, IR64 for O. sativa indica, TOG5681 for O. glaberrima, B88 for O. barthhi and W1654 for O. rufipogon) and were grown in June in the greenhouse at IRD, Montpellier, in Long Day condition (14h light/10h-dark during 6 to 8 weeks at 30°C, 80% relative humidity) and then transferred to Short Day condition (10h light/10h-dark at 30°C, 80% relative humidity) to induced flowering. Evaluation of panicle traits was investigated for 9 panicles per accessions under the same growth conditions and same time. 


### Tissue collection and RNA sequencing
For expression analysis, around 15 immature panicles from each accession and from each stage were collected from 4 days to 15 days after floral induction in 3 biological replicates. For sample collection, leaves surrounding the young panicle were removed by hand and the reproductive tissue was cut with a sharp blade under a stereomicroscope to confirm developmental stage. 
Immature panicle stages were defined as: Stage 1, Rachis Meristem (RM), Stage 2 Branching Meristem BM (i.e. initiation of primary branches, panicles with elongated primary and higher order branch development), Stage 3 Spikelet differentiation (SpM) and Stage 4 Floret differentiation and early flower organ development (**Fig. Mat Veg**). 
The reproductive tissues were frozen in liquid nitrogen and total RNAs (including small RNAs) were extracted using RNeasy Plant Mini Kit with RLT and RWT buffers (Qiaegn, France). DNAse treatments were performed using the RNAeasy-free DNase set (Qiagen, France). The RNA integrity numbers (RINs) of the extracted RNA measured by Agilent 2100 Bioanalyzer varied from 8,6 o 10). Stage and meristem specificity were validated with quantitative RT-PCR using stage specific markers genes. Primer sequences are available in Supplemental Table XX.
400ng of total RNA from stages 2 and 3 (from the 3 biological replicates) were used for performing libraries using the TrueSeq Stranded Total RNA with Ribo-Zero Plant prep Kit (Illumina) and following the standard protocol.  Libraries of 180 pb were quantified using a Bioanlyzer 2100 (Agilent) and paired-end 125 base sequencing was performed on the Get platform (Toulouse) using HiSeq2500. Samples were randomized across sequencing flow cells in 5 lanes. 
Add a table with reads nb for each librairy, Average input read length/%unique mapped reads/ %reads mapped to multiple loci/ final mapping  and % of reads on ribosomal seq (not removed from the ribodepletion kit)?

### Read alignment and Gene expression quantification

### RT-qPCR
First-strand cDNA was synthesized from 1μg of DNase-treated total RNA using the SuperScriptIII cDNA First-stand synthesis system (Invitrogen) as recommended by the manufacturer and was stored at -20°C. Large-scale RT-qPCR was used to study the expression by taking advantage of the Microfluidic Dynamic Array and according to the manufacturer’s instructions  (Fluidigm, San Diego, CA, USA).
Before performing real-time PCR, the sample mixture and assay mixture were prepared individually. 96 × 96 Dynamic Array Integrated Fluidic Circuit was loaded with cDNAs (sample inlets) and primer combinations (assay inlets) after specific target amplification (15 cycles) (STA) and exonuclease I treatment. A fast cycling protocol and EvaGreen (Bio-Rad) as dye was used on a BioMark machine. The experiment was performed at the CIRAD platform (Montpellier, France) following the workflow provided by the manufacturer. Three biological replicates were performed for each sample. Data were normalized using 4 genes (LOC_Os06g11170, LOC_Os06g48970, LOC_Os01g16970, LOC_Os03g61680). 
Relative gene expression were determined using the….methods





