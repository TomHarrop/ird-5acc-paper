library(tidyverse)
library(pheatmap)
library(ggpubr) 

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

prepare_for_plot <- function(dats) {
  dats <- get_expression(locus_ids = dats$locus_id,
                            dds = dds)
  
  dats <- dats %>%
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
  
  rownames(dats) <- dats$locus_id
  return(dats)
}

top_pc5 <- prepare_for_plot(pc5[1:50, ])

pc5 <- pc5 %>% arrange(desc(PC5))
bot_pc5 <- prepare_for_plot(pc5[1:50, ])


tp <- pheatmap(top_pc5[, 2:ncol(top_pc5)], 
               color = colorRampPalette(c( "white", "blue4"))(50),
               cluster_cols = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9)

bt <- pheatmap(bot_pc5[, 2:ncol(bot_pc5)], 
               color = colorRampPalette(c( "white", "blue4"))(50),
               cluster_cols = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9)


svg("../fig/fig-genes-of-pc5.svg",
    width = 12, height = 8)
ggarrange(plotlist = list(tp[[4]],
                          bt[[4]]))
dev.off()
