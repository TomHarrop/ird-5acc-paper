library(tidyverse)
library(viridis)
library(pheatmap)
library(grid)
library(gridExtra) 
library(magrittr)


source("../analyze-fluidigm/helper_functions_fluidigm.R")  


# Load file ---------------------------------------------------------------



load(file = "../data/fluidigm-confim-sampling.Rdata")

qpcr_exp <- bind_rows(chip3_exp, chip4_exp) %>%
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj",
                                     "Ob", "Og",
                                     "Or", "Osi")))


# Impute 0/ND -------------------------------------------------------------

# What is a the minimal expression besides 0?
# btw 0 must be intended as Non Detected.

# set it to one order of magnitude smaller than
# the minimum that was analyzed for a given amplicon

impute_zero <- function(expr) {
  min_expr <- 
    expr[expr > 0] %>%
    min()
  return(min_expr/10)
} 


qpcr_imputed <- 
  qpcr_exp %>% 
  group_by(target_name) %>%
  mutate(expression = case_when(expression == 0 ~ impute_zero(expression),
                                TRUE ~ expression)) 
  

# Plot --------------------------------------------------------------------

p <- 
  qpcr_imputed %>%
  {lineplot_fluidigm(nm = "Sampling Check", dat = .)} +
  scale_y_log10() +
  labs(title = "Marker genes behave as expected [semi-log plot]")

pdf("../fig/suppl-fig-qpcr-confirms-sampling.pdf",
    height = 6, width = 6,
    paper = "a4")
p %>% print()
dev.off()


# OLD ---------------------------------------------------------------------

 
# pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling-free-y.pdf",
#     height = 9,
#     width = 9)
# lineplot_fluidigm(nm = "Sampling Check",
#                   dat = qpcr_imputed,
#                   expr_val = expression) +
#   scale_y_log10() +
#   facet_wrap(facets = c('target_name', "species"),
#              scales = "free_y") +
#   labs(title = "Marker genes behave as expected [axis y log10]")
# dev.off()
# 
# pdf("../fig/suppl-fig-04b-qpcr-confirms-sampling-scaled.pdf",
#     height = 6,
#     width = 6)
# lineplot_fluidigm(nm = "Sampling Check",
#                   dat = qpcr_exp,
#                   expr_val = expr_scale_gene_spec_0_10) +
#   scale_y_log10() +
#   labs(title = "Marker genes behave as expected [axis y log10]")
# dev.off()
