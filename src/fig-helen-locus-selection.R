load("../data/all-sep-deseq.Rdata")
library(pheatmap)
library(tidyverse)
source("helper-functions.R")

# Prepare Data ------------------------------------------------------------


dds <- readRDS("../data-raw/dds.Rds")

load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")

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

pcx <- pc_spc$x %>% 
  as.data.frame(.) %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam) %>%
  left_join(mapman) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))

helen_loci <- read_delim(file = "../data-raw/LocusCluster3 to plot.csv",
                         delim = ";") %>%
  dplyr::rename(locus_id = "MSU_ID")

# Plot selected -----------------------------------------------------------

pdf("../fig/fig-helen-locus-selection.pdf",
    height = 30, 
    width = 20)
p <- helen_loci$locus_id %>%
  # as_factor(.) %>%
  get_expression(dds) %>%
  left_join(tf_fam)%>%
  left_join(annos) %>%
  left_join(mapman) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id",
                        "symbol",
                        "Family",
                        "DESCRIPTION"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 50,
                                       multi_line = T))
print(p)
dev.off()
4
