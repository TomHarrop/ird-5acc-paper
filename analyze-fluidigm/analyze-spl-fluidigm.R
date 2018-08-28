library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")

fluidigms <- c(spl_exp = "../data-raw/fluidigm_ok/SPok_CHIP4.xlsx",
               ref_exp = "../data-raw/fluidigm_ok/HK_CHIP4.xlsx") %>%
  map(., read_fluidigm)


norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

spl_exp <- fluidigms$spl_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

spl_heat <- spl_exp %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()


# Plot heatmap ------------------------------------------------------------

pdf("../fig/spl_fluidigm.pdf")
pheatmap(mat = spl_heat %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = spl_heat$locus_id)
dev.off()

