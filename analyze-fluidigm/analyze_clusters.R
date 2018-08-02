library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")


# load fluidigms -------------------------------------------------------

fluidigms <- c(clust_exp = "../data-raw/fluidigm/VariousCluster_CHIP3.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# Merge results -----------------------------------------------------------

chip3_exp <- fluidigms %>% select(-norms)

# merge normalizers and estimate expression -------------------------------

clust_exp <- fluidigms$clust_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

clust_heat <- clust_exp %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()


# Plot heatmap ------------------------------------------------------------

pdf("../fig/cluster_fluidigm.pdf",
    width = 12)
pheatmap(mat = clust_heat %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = clust_heat$locus_id)
dev.off()


# Plot both RNAseq and fluidigm -------------------------------------------

library(DESeq2)
library(gridExtra)
library(grid)

source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

clust_rnaseq <- unique(clust_exp$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "rufipogon", 
                                     "glaberrima", 
                                     "indica"))) 

clust_exp <- clust_exp %>% 
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj", 
                                     "Ob",
                                     "Or", 
                                     "Og", 
                                     "Osi"))) 


  
# id_fluidigm <- "LOC_Os06g47150"
# rnaseq_dat <- clust_rnaseq

# plot_both("LOC_Os04g36054",
#           fluidigm_dat = clust_exp,
#           rnaseq_dat = clust_rnaseq)

pdf(file = "../fig/cluster-fluidigm-rnaseq.pdf",
    height = 7)
unique(clust_exp$locus_id) %>% map(~plot_both(.,
                                             fluidigm_dat = clust_exp,
                                             rnaseq_dat = clust_rnaseq))
dev.off()
