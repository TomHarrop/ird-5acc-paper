# Note

**Do not base your work on commits in the `devel` branches, because we might amend them**

# Reproduce the analysis

To generate the analysis, you have to place these files into the `data-raw/` folder :

```
data-raw/
├── AP2_data9696n°2.xlsx
├── dds.Rds
├── famInfo.table.txt
├── geneInfo.table.txt
├── geneKeyword.table.txt
├── ListAP2_OKFLUIDGM.xlsx
├── MAPMAN BIN-Osa_MSU_v7.xlsx
├── OsOgObOrPTRAPdata_PaperTom.txt
└── Phenotype_PanicleSequenced_corrected.xlsx
```

Moreover you need to install these R packages from CRAN:

- `tidyverse`,
- `ggfortify`,
- `ggrepel`
- `readxl`,
- `gridExtra`,
- `leaps`
- `pheatmap`

And these from Bioconductor:

- `DESeq2`,
- `fgsea`

You might have to initiate some folder manually. This is my folder structure for the compiled paper.

```
ird-5acc-paper
├── bib
├── css
├── data
├── data-raw
├── docx
├── fig
├── html
├── ms
├── notes
├── src
├── tables
└── template
```
