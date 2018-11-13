library(tidyverse)
library(readxl)
library(DESeq2)

source("helper_functions_fluidigm.R")

fluidigms <- c(hb_exp = "extra-data/Hb_Chip4_domestication.xlsx",
               ref_exp = "extra-data/HK_CHIP4.xlsx") %>%
  map(., read_fluidigm)

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# merge normalizers and estimate expression -------------------------------

hb_fluid <- fluidigms$hb_exp %>%
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

hb_rnaseq <-
  unique(hb_fluid$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "glaberrima", 
                                     "rufipogon", 
                                     "indica"))) 

hb_fluid <- 
  hb_fluid %>% 
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

pdf(file = "../fig/hb-fluid-rnaseq-new.pdf",
    height = 7)
hb_fluid$locus_id %>% 
unique() %>% 
  map(~plot_both(.,
                 fluidigm_dat = hb_fluid,
                 rnaseq_dat = hb_rnaseq))
dev.off()
