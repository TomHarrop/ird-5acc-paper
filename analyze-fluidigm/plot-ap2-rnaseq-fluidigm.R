library(tidyverse)
library(readxl)
library(DESeq2)

source("helper_functions_fluidigm.R")
source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

# load fluidigms -------------------------------------------------------

fluidigms <- c(ap2_exp = "../data-raw/fluidigm/AP2_CHIP3.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)


# Analyze normalizers -----------------------------------------------------

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


# load rnaseq -------------------------------------------------------------

ap2_rnaseq <- unique(ap2_fluid$locus_id) %>%
  get_expression(dds = dds)
