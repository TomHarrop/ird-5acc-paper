load("../data/all-sep-deseq.Rdata")
library(pheatmap)
library(tidyverse)
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
  mutate(Family = ifelse(is.na(Family), "none", Family))

pdf("../fig/fig_ap2_from_pc1.pdf",
    width = 20,
    height = 60)
pcx %>%
  arrange(PC1) %>%
  filter(Family == "AP2-EREBP") %>%
  .$locus_id %>%
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

dev.off()

pdf("../fig/fig_mads_from_pc1.pdf",
    width = 20,
    height = 40)
pcx %>%
  arrange(PC1) %>%
  filter(Family == "MADS") %>%
  .$locus_id %>%
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

dev.off()

pdf("../fig/fig_zf_hd_from_pc1.pdf",
    width = 20,
    height = 25)
pcx %>%
  arrange(PC1) %>%
  filter(Family == "zf-HD") %>%
  .$locus_id %>%
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

dev.off()

table(tf_fam$Family)
