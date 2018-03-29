To generate the analysis, you have to place these files into the `data/` folder :

- `dds.Rds`,
- `MAPMAN BIN-Osa_MSU_v7.xlsx`,
- `Phenotype_PanicleSequenced_corrected..xlsx` (phenotypes, measured in Montpellier),
- `OsOgObOrPTRAPdata_PaperTom.txt` (phenotypes, measured in Cali)

Moreover you need to install these R packages from CRAN:

- `tidyverse`,
- `ggfortify`,
- `readxl`,
- `gridExtra`,
- `leaps`

And these from Bioconductor:

- `DESeq2`,
- `fgsea`
