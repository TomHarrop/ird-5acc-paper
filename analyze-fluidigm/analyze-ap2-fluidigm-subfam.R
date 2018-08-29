library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")

fluidigms <- c(ap2_exp = "../data-raw/fluidigm_ok/AP2ok_CHIP3.xlsx",
               ref_exp = "../data-raw/fluidigm_ok/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)


norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

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
  prepare_for_heat(scaled_val = expr_scale_gene_spec) %>%
  filter(complete.cases(.))



# Plot heatmap ------------------------------------------------------------

# pdf("../fig/hb_fluidigm.pdf")
# pheatmap(mat = hb_heat %>% select(-locus_id),
#          # scale = "row",
#          # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
#          # color = colorRampPalette(c( "white", "blue4"))(50),
#          color = viridis_pal()(50),
#          cluster_cols = F,
#          gaps_col = 1:4*5,
#          cellwidth = 9,
#          cellheight = 9,
#          labels_row = hb_heat$locus_id)
# dev.off()

# read subfamilies --------------------------------------------------------

ap2_sharoni_path <- "../data-raw/ap2s-sharoni/Supplementary_Table_S1.xls"

ap2_sharoni <- bind_cols(read_excel(ap2_sharoni_path,
                                    range = "A4:A166",
                                    col_names = "locus_id"),
                         read_excel(ap2_sharoni_path,
                                    range = "N4:O166",
                                    col_names = c("subgroub", "subfamily"))) %>%
  mutate(locus_id2 = paste0("LOC_", locus_id),
         subfam = paste(subfamily, subgroub, sep = " - ")) %>%
  select(locus_id2, subfamily) %>%
  mutate(subfam = subfamily) %>%
  select(-subfamily)

# %>%
#   select(class, locus_id2)


# heatmap with subfamilies ------------------------------------------------

ap2_heat_sfam <- ap2_heat %>%
  mutate(locus_id2 = str_split_fixed(string = locus_id,
                                     pattern = " ",
                                     2)[, 1]) %>%
  inner_join(ap2_sharoni) %>%
  mutate(subfam = as_factor(subfam)) %>%
  # filter(complete.cases(.)) %>%
  column_to_rownames("locus_id2")

pdf("../fig/ap2_fluidigm_families.pdf", width = 12)
pheatmap(mat = ap2_heat_sfam %>% select(-locus_id, -subfam),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat_sfam$locus_id,
         annotation_row = ap2_heat_sfam[, "subfam", drop = FALSE])
dev.off()


# line plot ap2 -----------------------------------------------------------

pdf("../fig/ap2-fludigm-lineplot.pdf",
    height = 50,
    width = 10)
ap2_exp %>%
  scale_tidy_fluidigm() %>%
  # pull(stage) %>%
  # as_factor() %>% 
  # as.numeric()
  lineplot_fluidigm()
dev.off()

