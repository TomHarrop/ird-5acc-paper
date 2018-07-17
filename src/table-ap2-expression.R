library(tidyverse)
library(DESeq2)
source("helper-functions.R")

# Read files --------------------------------------------------------------

dds <- readRDS(file = "../data-raw/dds.Rds")

ap2s <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id)) %>%
  filter(Family == "AP2-EREBP")


# Select Ap2s -------------------------------------------------------------

ap2_expr <- ap2s$locus_id %>%
  get_expression(dds) %>%
  select(locus_id, sample, `normalized expression`) %>%
  spread(key = sample,
         value = `normalized expression`) 

ap2_expr %>% write_csv("../tables/ap2-expression.csv")
