# load packages and data --------------------------------------------------

library(tidyverse)
library(viridis)
library(DESeq2)
library(pheatmap)


source("../src/helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")


# Dreb from My's phylo tree -----------------------------------------------


rav <- c("LOC_Os01g49830",
         "LOC_Os01g04750",
         "LOC_Os01g04800")


# Plot --------------------------------------------------------------------

pdf("../fig/rav_expr.pdf",
    height = 4, width = 10)
rav %>%
  get_expression(dds) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr()
dev.off()
