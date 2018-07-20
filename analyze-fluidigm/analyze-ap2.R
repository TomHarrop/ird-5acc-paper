library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")


# load fluidigms -------------------------------------------------------

fluidigms <- c(ap2_exp = "../data-raw/fluidigm/AP2_CHIP3.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# merge normalizers and estimate expression -------------------------------

ap2_exp <- fluidigms$ap2_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

ap2_heat <- ap2_exp %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

# Plot heatmap ------------------------------------------------------------

pdf("../fig/ap2_fluidigm.pdf",
    height = 10)
pheatmap(mat = ap2_heat %>% select(-target_name),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:5*4,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat$target_name)

dev.off()

