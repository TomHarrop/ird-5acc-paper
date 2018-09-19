# Reproduce the analysis - requirements

## Raw Data

To generate the analysis you need to download the `data-raw` folder from Dropbox.


## R packages

Moreover you need to install these R packages from CRAN:

- `tidyverse`,
- `ggfortify`,
- `ggrepel`
- `readxl`,
- `gridExtra`,
- `leaps`
- `pheatmap`
- `ggpubr`

These from Bioconductor:

- `DESeq2`,
- `fgsea`

And these from GitHub with devtools:

- `oryzr` --- type `devtools::install_github("TomHarrop/oryzr")`

## Folder structure

You might have to initiate some folder manually. This is my folder structure for the compiled paper.

```
├── analyze-data
├── analyze-fluidigm
├── analyze-promoters
├── bib
├── css
├── data
├── data-raw
├── docx
├── fig
├── html
├── Makefile
├── ms
├── notes
├── old
├── old-data
├── old-data-raw
├── readme.md
├── seq
├── src
├── tables
└── template
```

# Reproduce the analysis - transcriptomes

Run `make analysis` to analyze data, `make figures` to reproduce the figures and `make manuscript` to reproduce the mansucript.

Otherwise run the single individual scripts in folders.

# Parameters

Good genes (conserved) :

PC5 > .003
PC5 < .003

Good genes (DESeq2)

padj < 0.1 ?

# Graphical parameters

This are the graphical parameters that I am trying to use right now.

### Categorical variables

```
color_palette <- c("blue", "goldenrod")
```

### Continuous variables and scales

```
library(viridis)
viridis_pal()(50)
```

### Order of species

Reference - African (wild - domesticated) - Asian (wild - domesticated)

```
species <- factor(species,
                  levels = c("O. rufipogon",
                             "O. indica",
                             "O. japonica",
                             "O. barthii",
                             "O. glaberrima"))
```



# Note

**Do not base your work on commits in the `devel` branches, because we might amend them**
