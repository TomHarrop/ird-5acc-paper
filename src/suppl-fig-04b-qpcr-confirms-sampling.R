library(tidyverse)
# library(viridis)
# library(pheatmap)
# library(grid)
# library(gridExtra) 
library(magrittr)
library(cowplot)


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
  filter(target_name != "OsMADS14") %>% 
  {lineplot_fluidigm(nm = "Sampling Check", dat = .)} +
  scale_y_log10() +
  labs(title = "Marker genes behave as expected [semi-log plot]")


# Add tiff figure to plot -------------------------------------------------

# Note
#
# compressed with:
# tiffcp -c zip \
#     data-raw/FigMeristemCollect.tif \
#     data-raw/FigMeristemCollect-compr.tif 

p_img <- ggdraw() +  
  draw_image("../data-raw/FigMeristemCollect-compr.tif")

p_comb <- plot_grid(p_img, p,
                    labels = c("1.", "2."),
                    nrow = 2,
                    rel_heights = c(1, 2)) %>%
  add_sub(str_wrap("1. Developmental stages of immature panicles
                   collected for expression analysis. 
                   Stage 1: rachis meristem;
                   Stage2: formation of primary branch meristems,
                   elongation of primary branch meristem and formation
                   of axillary meristem; 
                   Stage3, spikelet meristem and floret differentitaion;
                   Stage 4, floral organ differenciation/development.
                   M, axillary merisyem, Fl, flower; Sp, spikelet,
                   RM, Rachis meristem; PbM, primary branch meristem;
                   ePbM, primary branch elongated;
                   Flm, floret meristem; St, stamen;
                   p, palea; l, lemma. 
                   2. In our RNA-seq samples, five selected marker genes
                   are expressed as expected.
                   This plot represents the fluidigm qPCR expression values
                   for each RNAseq sample relative and samples from two additional
                   developmental stages (In the RNAseq we have sequenced the samples
                   from developmental stages 2 and 3).
                   Each expression value is relative to the geometric mean expression
                   of 4 normalizers."),
          size = 11)

pdf("../fig/suppl-fig-qpcr-confirms-sampling.pdf",
    height = 10, width = 6.2,
    paper = "a4")

p_comb %>%
  ggdraw() %>%
  print()

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
