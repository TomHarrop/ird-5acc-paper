library(tidyverse)
library(viridis)
library(pheatmap)
library(gridExtra)
source("../analyze-fluidigm/helper_functions_fluidigm.R")  


load("../data/ap2_fluidigm.R")
load("../data/hb_fluidigm.R")

cls <- read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  dplyr::rename(locus_id = "MsuID")

cl4 <- cls %>% filter(cluster == 4)
cl5 <- cls %>% filter(cluster == 5)

ap2_heat5 <- ap2_exp %>%
  filter(locus_id %in% cl5$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c5a <- pheatmap(mat = ap2_heat5 %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50), main = "AP2 cluster5",
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat5$locus_id)


ap2_heat4 <- ap2_exp %>%
  filter(locus_id %in% cl4$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c4a <- pheatmap(mat = ap2_heat4 %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50), main = "AP2 cluster4",
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat4$locus_id)

hb_heat4 <- hb_exp %>%
  filter(locus_id %in% cl4$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c4h <- pheatmap(mat = hb_heat4 %>% select(-locus_id),
                # scale = "row",
                # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
                # color = colorRampPalette(c( "white", "blue4"))(50),
                color = viridis_pal()(50), main = "HB cluster4",
                cluster_cols = F,
                gaps_col = 1:4*5,
                cellwidth = 9,
                cellheight = 9,
                labels_row = hb_heat4$locus_id)

pdf("../fig/fig-06-fludigm-ap2-hb-DRAFT.pdf")
grid.arrange(c4h[[4]], c4a[[4]], c5a[[4]])
dev.off()
