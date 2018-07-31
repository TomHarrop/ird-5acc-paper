# load packages and data --------------------------------------------------

library(tidyverse)
library(viridis)
library(DESeq2)
library(pheatmap)


source("../src/helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

# Read AP2 ids-------------------------------------------------------------

ap2s <- read_delim("locus-id-ordered.txt",
           delim = " ",
           col_names = FALSE) %>%
  mutate(locus_id = paste(X1, X2, sep = "_")) %>%
  select(locus_id)


# Get expression ----------------------------------------------------------


expr <- ap2s %>%
  pull(locus_id) %>%
  get_expression(dds) %>%
  mutate(locus_id = as_factor(locus_id))

to_heat <- expr %>%
  group_by(locus_id, species, stage) %>%
  summarise(to_plot = median(by_locus_species)) %>%
  ungroup() %>%
  mutate(stage_species = paste(stage, species, sep = "_")) %>%
  select(locus_id, stage_species, to_plot) %>%
  spread(key = stage_species, value = to_plot) %>%
  as.data.frame() %>%
  mutate(locus_id_fact = as_factor(locus_id)) %>%
  column_to_rownames("locus_id") 

# Plot heatmap ------------------------------------------------------------

p <- pheatmap(to_heat %>% select(-locus_id_fact),
              color = viridis_pal()(50),
              cluster_cols = F,
              cluster_rows = F,
              gaps_col = 5,
              cellwidth = 12,
              cellheight = 9)

pdf("../fig/heatmap-ap2-order-phylogen-tree.pdf", height = 20)
print(p)
dev.off()
