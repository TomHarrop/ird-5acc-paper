library(tidyverse)
library(pheatmap)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

pc5 <- pcro %>%
  select(PC5, locus_id) %>%
  arrange(PC5)

top_pc5 <- pc5[1:50, ]

top_pc5 <- get_expression(locus_ids = top_pc5$locus_id,
                          dds = dds)

top_pc5 <- top_pc5 %>%
  group_by(locus_id, species, stage) %>%
  # mutate(scaled_expression = scale(`normalized expression`)) %>%
  summarise(by_locus_species = median(by_locus_species)) %>%
  ungroup() %>%
  mutate(stage_species = paste(stage, species, sep = "_")) %>%
  select(locus_id, stage_species, by_locus_species) %>%
  spread(key = stage_species, value = by_locus_species) %>%
  left_join(annos) %>%
  mutate(locus_id = paste(locus_id, symbol, sep = "_")) %>%
  select(-symbol) %>%
  as.data.frame(.)


rownames(top_pc5) <- top_pc5$locus_id
svg("../fig/fig-genes-of-pc5.svg",
    width = 5, height = 8)
pheatmap(top_pc5[, 2:ncol(top_pc5)], 
         color = colorRampPalette(c( "white", "blue4"))(50),
         cluster_cols = F,
         gaps_col = 5,
         cellwidth = 9,
         cellheight = 9)
dev.off()
