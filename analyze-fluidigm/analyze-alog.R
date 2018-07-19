library(tidyverse)
library(readxl)


# Function for Geometric Mean ---------------------------------------------

gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Function that reads fluidigm --------------------------------------------

read_fluidigm <- function(path) {
  out <- read_excel(path) %>%
    select(sampleName, GeneName,
           LOC_Name, Value, Type__1) %>%
    rename(sample_name = sampleName,
           target_name = GeneName,
           locus_id = LOC_Name,
           ct_value = Value,
           type = Type__1) %>%
    mutate(ct_value = as.numeric(ct_value))
  out 
}


# load fluidigms -------------------------------------------------------

fluidigms <- c(alog_exp = "../data-raw/ALOG_CHIP3.xlsx",
               ref_exp = "../data-raw/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# merge normalizers and estimate expression -------------------------------

alog_exp <- fluidigms$alog_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

# Store descriptive names for stages

stage_names <- c(N1 = "Rachis Meristem",
                 N2 = "Primary Branch Meristem",
                 N3 = "Spikelet meristem",
                 N4 = "Young Flower Development")

# Store descriptive names for species

species_names <- c(B = "Barthii",
                   G = "Glaberrima",
                   I = "Indica",
                   J = "Japonica",
                   R = "Rufipogon")

# Scale gene expression and scale gene expression separating 

alog_exp <- alog_exp %>%
  # make a variable that encodes the species
  mutate(species = str_sub(sample_name, 1, 1)) %>%
  # make a variable that encodes for the stage
  mutate(stage = str_sub(sample_name, 3, 4)) %>%
  # make a variable that encode replicate
  mutate(replicate = str_sub(sample_name, 6, 6)) %>%
  # scale expression by gene
  group_by(target_name) %>%
  mutate(expr_scale_gene = scale(expression)) %>%
  ungroup() %>%
  # scale expression by geme end species
  group_by(target_name, species) %>%
  mutate(expr_scale_gene_spec = scale(expression)) %>%
  ungroup() %>%
  mutate(expr_scale_gene_spec = ifelse(is.na(expr_scale_gene_spec), 0, expr_scale_gene_spec)) %>%
  # last, save descriptive names for stages and species
  mutate(stage = stage_names[stage]) %>%
  mutate(stage = factor(stage, 
                        levels = c("Rachis Meristem",
                                   "Primary Branch Meristem",
                                   "Spikelet meristem",
                                   "Young Flower Development"))) %>%
  mutate(species = species_names[species])

