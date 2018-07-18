library(tidyverse)
load("../data/rlog-pca.Rdata")

# Family def from funricegene
fams <- read_delim("../data-raw/famInfo.table.txt",
                   delim = "\t") %>%
  select(MSU, Name) %>%
  dplyr::rename(locus_id = "MSU",
                Family = "Name") %>%
  filter(locus_id != "none")

# merge with PC5
dats <- pcro %>% 
  left_join(fams) %>%
  arrange(desc(PC5)) %>%
  mutate(rank_pc = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family)) %>%
  select(locus_id, PC5, Family) %>%
  distinct() %>%
  # cutof of .003
  mutate(de = abs(PC5) > .003)

# estimate enrichement

# Define phyper wrapper that contains "..."
# So that it can be used in pmap with extra variables
phyper2 <- function(q, m, n, k, ...) phyper(q, m, n, k, lower.tail = FALSE)


# Test enrichment
# inspired from
# https://github.com/GuangchuangYu/DOSE/blob/master/R/enricher_internal.R
tst <- dats %>%
  # white balls
  add_count(Family) %>%
  dplyr::rename(m = "n") %>%
  # balls drawn
  add_count(de) %>%
  dplyr::rename(k = "n") %>%
  # black balls
  add_tally() %>%
  mutate(n = n - m) %>%
  # white balls drawn
  filter(de) %>%
  group_by(Family) %>%
  mutate(q = n()) %>%
  select(-locus_id, -PC5) %>%
  distinct() %>%
  ungroup() %>%
  # estimate pval
  mutate(pval = pmap(., .f = phyper2, lower.tail = FALSE)) %>%
  mutate(pval = as.numeric(pval)) %>%
  arrange(pval)

fams %>%
  filter(!locus_id %in% dats$locus_id)
