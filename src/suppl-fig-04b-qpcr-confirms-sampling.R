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
                  dat = qpcr_exp,
                  expr_val = expression) +
  scale_y_log10() +
  facet_wrap(facets = c('target_name', "species"),
             scales = "free_y")

pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling.pdf",
    height = 6, width = 6)
lineplot_fluidigm(nm = "Sampling Check",
                  dat = qpcr_exp) +
  scale_y_log10() +
  labs(title = "Marker genes behave as expected [axis y log10]")
dev.off()

pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling-free-y.pdf",
    height = 9,
    width = 9)
lineplot_fluidigm(nm = "Sampling Check",
                  dat = qpcr_exp,
                  expr_val = expression) +
  scale_y_log10() +
  facet_wrap(facets = c('target_name', "species"),
             scales = "free_y") +
  labs(title = "Marker genes behave as expected [axis y log10]")
dev.off()

pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling-scaled.pdf",
    height = 6,
    width = 6)
lineplot_fluidigm(nm = "Sampling Check",
                  dat = qpcr_exp,
                  expr_val = expr_scale_gene_spec_0_10) +
  scale_y_log10() +
  labs(title = "Marker genes behave as expected [axis y log10]")
dev.off()
