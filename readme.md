To generate the analysis, you have to place these files into the `data-raw/` folder :

- `dds.Rds`,
- `MAPMAN BIN-Osa_MSU_v7.xlsx`,
- `Phenotype_PanicleSequenced_corrected..xlsx` (phenotypes, measured in Montpellier),
- `OsOgObOrPTRAPdata_PaperTom.txt` (phenotypes, measured in Cali)
- `AP2_data9696n°2.xlsx` Fluidigm on selected AP2 genes
- `ListAP2_OKFLUIDGM.xlsx` in fluidigm, which gene gives best results?

Moreover you need to install these R packages from CRAN:

- `tidyverse`,
- `ggfortify`,
- `readxl`,
- `gridExtra`,
- `leaps`
- `pheatmap`

And these from Bioconductor:

- `DESeq2`,
- `fgsea`
