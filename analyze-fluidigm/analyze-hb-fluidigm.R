library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")

fluidigms <- c(hb_exp = "../data-raw/fluidigm_ok/HBok_CHIP4.xlsx",
               ref_exp = "../data-raw/fluidigm_ok/HK_CHIP4.xlsx") %>%
  map(., read_fluidigm)


norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

hb_exp <- fluidigms$hb_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

hb_heat <- hb_exp %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat() %>%
  filter(complete.cases(.))
  


# Plot heatmap ------------------------------------------------------------

pdf("../fig/hb_fluidigm.pdf")
pheatmap(mat = hb_heat %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = hb_heat$locus_id)
dev.off()

# read subfamilies --------------------------------------------------------

hb_sfam <- read_csv("../data-raw/hb_genes_jain2008.csv") %>%
  dplyr::rename(locus_id2 = "msuId") 

# %>%
#   select(class, locus_id2)
  

# heatmap with subfamilies ------------------------------------------------

hb_heat_sfam <- hb_heat %>%
  mutate(locus_id2 = str_split_fixed(string = locus_id,
                                   pattern = " ",
                                   2)[, 1]) %>%
  inner_join(hb_sfam) %>%
  mutate(class = as_factor(class),
         locus_id = paste(locus_id2, symbol)) %>%
  select(-symbol) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("locus_id2")

pdf("../fig/hb_fluidigm_families.pdf", width = 12)
pheatmap(mat = hb_heat_sfam %>% select(-locus_id, -class),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = hb_heat_sfam$locus_id,
         annotation_row = hb_heat_sfam[, "class", drop = FALSE])
dev.off()
