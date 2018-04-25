# Set up environment ------------------------------------------------------

library(tidyverse)
# library(kableExtra)
library(fgsea)
library(DESeq2)
# library(knitr)
library(gridExtra)
load("../data/rlog-pca.Rdata")
load("../data/mapman.Rdata")
source("helper-functions.R")


# estimate enrichment for PC5 ---------------------------------------------

mapman_list <- c(split(mapman$IDENTIFIER, mapman$level1),
                 split(mapman$IDENTIFIER, mapman$level2),
                 split(mapman$IDENTIFIER, mapman$level3),
                 split(mapman$IDENTIFIER, mapman$level4),
                 split(mapman$IDENTIFIER, mapman$level4),
                 split(mapman$IDENTIFIER, mapman$level6),
                 split(mapman$IDENTIFIER, mapman$level7))

gsea_pc5 <- test_gsea(set_names(pcro$PC5, nm = pcro$locus_id), mapman = mapman_list)

svg("../tables/table-pc5-mapman.svg",
    width = 10,
    height = 2)
gsea_pc5 %>% select(pathway, padj) %>%
  filter(pathway %in% c("RNA.regulation of transcription.AP2/EREBP, APETALA2/ethylene-responsive element binding protein family",
                        "RNA.regulation of transcription.C2C2(Zn) Constans-like zinc finger family (CO-like)",
                        "RNA.regulation of transcription.general transcription",
                        "RNA.regulation of transcription.MADS box transcription factor family",
                        "RNA.regulation of transcription.MYB domain transcription factor family")) %>%
  grid.table()
dev.off()
  # kable("html") %>%
  # kable_styling(bootstrap_options = c("striped", "hover")) %>%
  # # write("../tables/table-pc5-mapman.html")
