library(tidyverse)
library(fgsea)
library(magrittr)

dds <- readRDS("../data-raw/dds.Rds")
source("helper-functions.R")
filter_abs_pc <- .00185

set.seed(1)

# Load PCA and TF families ------------------------------------------------

load("../data/rlog-pca.Rdata")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

pcx <- pcro %>%
  as.data.frame(.) %>%
  left_join(tf_fam) %>%
  arrange(desc(PC5)) %>%
  mutate(rank_pc5 = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))


# Test Enrichment ---------------------------------------------------------

gsea_pc5 <- test_gsea(set_names(pcx$PC5, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
  dplyr::rename(Family = "pathway") %>%
  mutate(leadingEdge = leadingEdge %>% 
           map(~paste(., collapse = " - ")) %>%
           unlist()) %>%
  as_tibble() %T>%
  write_csv(path = "../tables/PC5-TF-enrichment-gsea.csv")



