library(tidyverse)

load("../data/all-sep-deseq.Rdata")
dds <- readRDS("../data-raw/dds.Rds")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

pcx <- pc_spc$x %>% 
  as.data.frame(.) %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam) %>%
  arrange(desc(PC1)) %>%
  mutate(rank_pc1 = 1:nrow(.)) 

pcx %>% top_n(15) %>% write_csv("../tables/pc1-last-15.csv")

pcx %>%
  arrange(PC1) %>%
  top_n(15) %>%
  write_csv("../tables/pc1-top-15,csv")
