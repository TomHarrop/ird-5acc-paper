
## Materials and Methods

> n.b. I haven't added any changes from this section to the MS file

### Plant Material and Growth Conditions

Panicle traits diversity were investigated in accessions (**Table S1 Plantinfo**;) of wild and domesticated, Asian and African rice.  O. sativa, O. glaberrima, O. barthii  and o. plants were grown in Cali (Colombia) and in Montpellier (France) for O. rufipogon.
At panicle maturity, the 3 main panicle from 3 plants per accession per repeat were collected (i.e. 18 panicles/ accession). Each panicle was spread out and fixed on a white paper by tape. Panicle structure and seed number were analyses using P-TRAP software (Al-Tam et al., 2013). The quantified panicle traits included rachis length (RL in cm), total panicle length (PanL in cm), primary branch number (PbN), primary branch average length (PbL in cm), primary branch internode average length (PbintL in cm), secondary branch number (SbN), secondary branch average length (SBL in cm), secondary branch internode average length (SbintL in cm) and spikelet number (SpN). Phenotype description and statistical analysis of variance were performed by using functions in R software.
For expression and sequencing analysis, 5 accessions have been used (Niponbarre for O. sativa japonica, IR64 for O. sativa indica, TOG5681 for O. glaberrima, B88 for O. barthhi and W1654 for O. rufipogon) and were grown in June in the greenhouse at IRD, Montpellier, in Long Day condition (14h light/10h-dark during 6 to 8 weeks at 30°C, 80% relative humidity) and then transferred to Short Day condition (10h light/10h-dark at 30°C, 80% relative humidity) to induced flowering. Evaluation of panicle traits was investigated for 9 panicles per accessions under the same growth conditions and same time.


### Tissue collection and RNA sequencing
For expression analysis, around 15 immature panicles from each accession and from each stage were collected from 4 days to 15 days after floral induction in 3 biological replicates. For sample collection, leaves surrounding the young panicle were removed by hand and the reproductive tissue was cut with a sharp blade under a stereomicroscope to confirm developmental stage.
Immature panicle stages were defined as: Stage 1, Rachis Meristem (RM), Stage 2 Branching Meristem BM (i.e. initiation of primary branches, panicles with elongated primary and higher order branch development), Stage 3 Spikelet differentiation (SpM) and Stage 4 Floret differentiation and early flower organ development (**Fig. Mat Veg**).

>here: From the exact same plants that were sampled for the RNAseq ?

The reproductive tissues were frozen in liquid nitrogen and total RNAs (including small RNAs) were extracted using RNeasy Plant Mini Kit with RLT and RWT buffers (Qiaegn, France). DNAse treatments were performed using the RNAeasy-free DNase set (Qiagen, France). The RNA integrity numbers (RINs) of the extracted RNA measured by Agilent 2100 Bioanalyzer varied from 8,6 o 10). Stage and meristem specificity were validated with quantitative RT-PCR using stage specific markers genes. Primer sequences are listed in Table S2 **supp-table-PrimerList**.
400ng of total RNA from stages 2 and 3 (from the 3 biological replicates) were used for performing libraries using the TrueSeq Stranded Total RNA with Ribo-Zero Plant prep Kit (Illumina) and following the standard protocol.  Libraries of 180 pb were quantified using a Bioanlyzer 2100 (Agilent) and paired-end 125 base sequencing was performed on the Get platform (Toulouse) using HiSeq2500. Samples were randomized across sequencing flow cells in 5 lanes.

### Read alignment and Gene expression quantification

### Gene set enrichment

GSEA: Subramanian 2005, https://doi.org/10.1073/pnas.0506580102

FGSEA: Sergushichev 2016, https://doi.org/10.1101/060012

The pvalue of GSEA are estimated empirically by permutation.

### RT-qPCR
First-strand cDNA was synthesized from 1μg of DNase-treated total RNA using the SuperScriptIII cDNA First-stand synthesis system (Invitrogen) as recommended by the manufacturer and was stored at -20°C. Large-scale RT-qPCR was used to study the expression by taking advantage of the Microfluidic Dynamic Array and according to the manufacturer’s instructions  (Fluidigm, San Diego, CA, USA).
Before performing real-time PCR, the sample mixture and assay mixture were prepared individually. 96 × 96 Dynamic Array Integrated Fluidic Circuit was loaded with cDNAs (sample inlets) and primer combinations (assay inlets) after specific target amplification (15 cycles) (STA) and exonuclease I treatment. A fast cycling protocol and EvaGreen (Bio-Rad) as dye was used on a BioMark machine. The experiment was performed at the CIRAD platform (Montpellier, France) following the workflow provided by the manufacturer. Three biological replicates were performed for each sample. Data were normalized using 4 genes (LOC_Os06g11170, LOC_Os06g48970, LOC_Os01g16970, LOC_Os03g61680).
Relative gene expression were determined using the….methods
