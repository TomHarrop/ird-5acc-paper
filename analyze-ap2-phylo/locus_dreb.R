# load packages and data --------------------------------------------------

library(tidyverse)
library(viridis)
library(DESeq2)
library(pheatmap)


source("../src/helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")


# Dreb from My's phylo tree -----------------------------------------------


dreb <- c("LOC_Os02g55380", "LOC_Os06g08340",
          "LOC_Os06g40150", "LOC_Os01g07120",
          "LOC_Os05g27930", "LOC_Os02g42585",
          "LOC_Os04g44670", "LOC_Os02g51670",
          "LOC_Os06g11860", "LOC_Os03g09170",
          "LOC_Os08g31580", "LOC_Os09g20350",
          "LOC_Os09g35010", "LOC_Os09g35020",
          "LOC_Os09g35030", "LOC_Os06g06970",
          "LOC_Os02g45450", "LOC_Os04g48350",
          "LOC_Os04g55520", "LOC_Os06g07030",
          "LOC_Os03g15660", "LOC_Os04g36640",
          "LOC_Os02g45420", "LOC_Os02g13710",
          "LOC_Os06g36000", "LOC_Os02g43940",
          "LOC_Os04g46400", "LOC_Os10g41130",
          "LOC_Os02g43970", "LOC_Os04g46440")



# Plot --------------------------------------------------------------------

pdf("../fig/dreb_expr.pdf",
    height = 16, width = 10)
dreb %>%
  get_expression(dds) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr()
dev.off()


# DREB nice clade ---------------------------------------------------------

dreb_clade <- c("LOC_Os02g13710",
                "LOC_Os06g36000",
                "LOC_Os02g43940",
                "LOC_Os04g46400",
                "LOC_Os10g41130",
                "LOC_Os02g43970",
                "LOC_Os04g46440")

pdf("../fig/dreb_clade_expr.pdf",
    height = 8, width = 10)
dreb_clade %>%
  get_expression(dds) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr()
dev.off()
