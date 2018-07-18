library(tidyverse)
library(readxl)
library(pheatmap)
source("helper_functions_fluidigm.R")

# load data
# from Tom's edi-hk-check/calculate_relq.R
headers <- c("Chamber", "Sample_Name", "Sample_Type", "Sample_rConc",
             "target_name", "target_type", "Ct_Value", "Ct_Calibrated_rConc",
             "Ct_Quality", "Ct_Call", "Ct_Threshold",
             "tm_in_range", "tm_out_range", "tm_peak_ratio")
dat <- read_excel(path = "../data-raw/AP2_data9696n°2.xlsx",
                  range = "A12:N8172") 
colnames(dat) <- headers

# How many measures per target?
table(dat$target_name)

# remove negative control (water?)
# No, some have Ct, what to do with them?
table(dat$Sample_Name)
summary(dat %>%
          filter(Sample_Name == "H20") %>%
          select(Ct_Value))
# for now remove.
dat <- dat %>% filter(Sample_Name != "H20")

# Remove standard curves?
# Do all sampl in standard curves contain the word "Mix" in their names? Yes
dat <- dat %>% filter(!grepl("Mix", Sample_Name))

# How many measurement now?
table(dat$target_name)
table(dat$Sample_Name)
table(dat$Sample_Type)


# Take geometric mean of normalizers ------------------------------------------------

# there are 3 normalizer per sample,
# you must take the geometric mean
norms <- c("ACT2 : Loc_Os11g06390",
           "HK09",
           "HK04")

# make a separate dataset for those
norms <- dat %>%
  # filter rows for normalizers
  filter(target_name %in% norms)

# how many measurements?
table(norms$target_name)
table(norms$Sample_Name)

# Take the geomoetric mean for every sample
# Check that the function for geometric mean works properly!
norms <- norms %>%
  group_by(Sample_Name) %>%
  summarize(norm_geom_mean = gm_mean(Ct_Value))



# Correct for normalizers and take exponential -------------------------

dat <- dat %>%
  # merge back normalizers to dataset
  left_join(norms, by = "Sample_Name") %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(Ct_Value - norm_geom_mean))) %>%
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

dat <- dat %>%
  # make a variable that encodes the species
  mutate(species = str_sub(Sample_Name, 1, 1)) %>%
  # make a variable that encodes for the stage
  mutate(stage = str_sub(Sample_Name, 3, 4)) %>%
  # make a variable that encode replicate
  mutate(replicate = str_sub(Sample_Name, 6, 6)) %>%
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

ap2_fluidigm <- dat
save(ap2_fluidigm, file = "../data/AP2-fluidigm.Rdata")


# Select only high quality genes ------------------------------------------

ap2_good <- read_excel("../data-raw/ListAP2_OKFLUIDGM.xlsx",
                       range = "A3:F55")

colnames(ap2_good) <- c("msu", "rap", "symbol", "primer_f", "primer_r", "ha_primers")

dat <- dat %>%
  filter(target_name %in% ap2_good$symbol)

unique(dat$target_name)


# reshape and plot --------------------------------------------------------

dat <- dat %>%
  group_by(target_name, species, stage) %>%
  summarise(median_expr = median(expr_scale_gene_spec)) %>%
  ungroup()

dat <- dat %>% 
  mutate(species_stage = paste(species, stage, sep = "_")) %>%
  select(target_name, species_stage, median_expr) %>% 
  spread(key = species_stage, value = median_expr) %>%
  as.data.frame(.)

rownames(dat) <- dat$target_name

# pheatmap(mat = dat[, 2:ncol(dat)],
#          # scale = "row",
#          # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
#          color = colorRampPalette(c( "white", "blue4"))(50),
#          cluster_cols = F,
#          gaps_col = 1:5*4,
#          cellwidth = 9,
#          cellheight = 9, 
#          filename = "../fig/figure-heatmap-fluidigm-ap2.pdf")
