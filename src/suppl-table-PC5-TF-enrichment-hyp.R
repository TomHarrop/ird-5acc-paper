library(tidyverse)
library(magrittr)
# library(fgsea)
dds <- readRDS("../data-raw/dds.Rds")
# source("helper-functions.R")
filter_abs_pc <- .00185

# Load PCA and TF families ------------------------------------------------

load("../data/rlog-pca.Rdata")

tf_fam <- 
  readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

# Prepare data - merge PC and TF ------------------------------------------

pcx <- 
  pcro %>%
  as.data.frame(.) %>%
  left_join(tf_fam) %>%
  arrange(desc(PC5)) %>%
  mutate(rank_pc5 = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))


# Test Enrichment ---------------------------------------------------------

genes <- 
  pcx %>% 
  select(PC5, locus_id, Family) %>%
  mutate(Family = na_if(Family, "none"),
         relevant = case_when(PC5 > filter_abs_pc ~ "yes",
                              PC5 < -filter_abs_pc ~ "yes",
                              TRUE ~ "no")) %>%
  filter(complete.cases(.))

make_contingency <- function(fam, data = genes) {
  tibble(family = fam,
         yes = data %>%
           filter(Family == fam, 
                  relevant == "yes") %>%
           nrow(),
         no = data %>%
           filter(Family == fam, 
                  relevant == "no") %>%
           nrow(),
         all_genes = genes %>%
           nrow,
         rels  = genes %>%
           filter(relevant == "yes") %>%
           nrow())
}

genelist <- 
  genes$Family %>%
  unique() %>%
  map(make_contingency) %>%
  reduce(bind_rows)

# Define phyper wrapper that contains "..."
# So that it can be used in pmap with extra variables
phyper2 <- function(q, m, n, k, ...) phyper(q, m, n, k, lower.tail = FALSE)

get_ids <- function(fam) {
  genes %>%
    filter(Family == fam) %$%
    locus_id %>%
    paste0(collapse = " - ")
}

# Test enrichment
# inspired from
# # https://github.com/GuangchuangYu/DOSE/blob/master/R/enricher_internal.R

cont <- 
  genelist %>%
  mutate(q = yes,  # white balls drawn
         m = no + yes, # White balls 
         n = all_genes - no - yes, # Black balls
         k = rels) %>% # balls drawn
  mutate(pval = pmap(., .f = phyper2, lower.tail = FALSE),
         pval = as.numeric(pval),
         padj = pval %>% p.adjust(method = "fdr"),
         genes = family %>% map(get_ids))  %>%
  arrange(pval) 
  
# Save
cont %>%
  mutate(genes = flatten_chr(genes)) %>%
  write_csv("../tables/tfs-of-pc5-hypergeom.csv")
