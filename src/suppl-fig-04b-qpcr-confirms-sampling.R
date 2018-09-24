library(tidyverse)
library(viridis)
library(pheatmap)
library(grid)
library(gridExtra)


source("../analyze-fluidigm/helper_functions_fluidigm.R")  


load(file = "../data/fluidigm-confim-sampling.Rdata")

qpcr_exp <- bind_rows(chip3_exp, chip4_exp) %>%
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj",
                                     "Ob", "Og",
                                     "Or", "Osi")))

lineplot_fluidigm(nm = "Sampling Check",
                  dat = qpcr_exp)

pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling.pdf",
    height = 6, width = 6)
lineplot_fluidigm(nm = "Sampling Check",
                  dat = qpcr_exp) +
  labs(title = "Marker genes behave as expected - more or less ")
dev.off()
