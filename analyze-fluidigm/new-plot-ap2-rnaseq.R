library(tidyverse)
library(readxl)
library(DESeq2)

source("helper_functions_fluidigm.R")

fluidigms <- c(ap2_exp = "extra-data/AP2_Fluidgm_Domestication list.xlsx",
               ref_exp = "extra-data/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# merge normalizers and estimate expression -------------------------------

ap2_fluid <- fluidigms$ap2_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Plot both RNAseq and fluidigm -------------------------------------------

library(DESeq2)
library(gridExtra)
library(grid)

source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

ap2_rnaseq <-
  unique(ap2_fluid$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "glaberrima", 
                                     "rufipogon", 
                                     "indica"))) 

ap2_fluid <- 
  ap2_fluid %>% 
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj", 
                                     "Ob",
                                     "Og", 
                                     "Or", 
                                     "Osi"))) 



# id_fluidigm <- "LOC_Os06g47150"
# rnaseq_dat <- clust_rnaseq

# plot_both("LOC_Os04g36054",
#           fluidigm_dat = clust_exp,
#           rnaseq_dat = clust_rnaseq)

pdf(file = "../fig/ap2-fluid-rnaseq-new.pdf",
    height = 7)
unique(ap2_fluid$locus_id) %>% 
  map(~plot_both(.,
                 fluidigm_dat = ap2_fluid,
                 rnaseq_dat = ap2_rnaseq))
dev.off()
