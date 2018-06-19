# Plot result of clustering


# Load Data ---------------------------------------------------------------


library(tidyverse)
library(readr)
dds <- readRDS("../data-raw/dds.Rds")
source("helper-functions.R")
load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
library(pheatmap)

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

mapman <- mapman %>%
  dplyr::rename(locus_id = "IDENTIFIER") %>%
  filter(!duplicated(locus_id))


clusts <- read_csv(file = "../data-raw/annotated_clusters_scaled_l2fc.csv") %>% 
  split(.$cluster)

plot_heatmap <- function(ids,
                         norm = by_locus_species,
                         cutree_rows = 1,
                         filter_abs_pc = .003) {
  norm <- enquo(norm)
  to_heat <- ids %>%
    get_expression(dds) %>%
    group_by(locus_id, species, stage) %>%
    summarise(to_plot = median(!!norm)) %>%
    ungroup() %>%
    mutate(stage_species = paste(stage, species, sep = "_")) %>%
    select(locus_id, stage_species, to_plot) %>%
    spread(key = stage_species, value = to_plot) %>%
    as.data.frame() %>%
    column_to_rownames("locus_id") 
  
  
  p <- pheatmap(to_heat,
                color = colorRampPalette(c( "white", "blue4"))(50),
                show_rownames = F,
                # fontsize = 5,
                cutree_cols = 2,
                cluster_cols = F,
                # cluster_rows = F,
                gaps_col = 5,
                # cutree_rows = cutree_rows,
                cellwidth = 9,
                cellheight = 5)
  
  return(p)
}

clusts %>% map(~plot_heatmap(.$MsuID))
