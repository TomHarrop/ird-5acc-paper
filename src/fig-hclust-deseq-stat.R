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

tst <- spc_res_df %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam)%>%
  left_join(annos) %>%
  left_join(mapman) %>%
  filter(Category == "TF") %>%
  mutate(annos = paste(locus_id, Family, symbol, DESCRIPTION, sep = "--"))


pdf("../fig/fig-hclust-deseq-stat.pdf", height = 210, width = 20)
pheatmap(tst %>%
           select(stat_japonica,
                  stat_barthii,
                  stat_glaberrima,
                  stat_rufipogon,
                  stat_indica),
         scale = "column", cellwidth = 10, cellheight = 10,
         treeheight_row = 200,
         labels_row = tst$annos,
         cluster_cols = F, cutree_rows = 20)
dev.off()

