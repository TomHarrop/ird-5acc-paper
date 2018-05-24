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
  mutate(Family = ifelse(is.na(Family), "none", Family))

# Plot selected -----------------------------------------------------------

plot_kword <- function(kword = "AP2-EREBP") {
  p <- pcx %>%
    arrange(PC1) %>%
    filter(Family == kword) %>%
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
  return(p)
}


pdf("../fig/fig_ap2_from_pc1.pdf",
    width = 20,
    height = 60)
print(plot_kword())
dev.off()

pdf("../fig/fig_mads_from_pc1.pdf",
    width = 20,
    height = 40)
print(plot_kword("MADS"))
dev.off()

pdf("../fig/fig_zf_hd_from_pc1.pdf",
    width = 20,
    height = 25)
print(plot_kword("zf-HD"))
dev.off()


pdf("../fig/fig_NAC_from_pc1.pdf",
    width = 20,
    height = 30)
print(plot_kword("NAC"))
dev.off()

pdf("../fig/fig_HB_from_pc1.pdf",
    width = 20,
    height = 60)
print(plot_kword("HB"))
dev.off()

pdf("../fig/fig_SPB_from_pc1.pdf",
    width = 20,
    height = 20)
print(plot_kword("SBP"))
dev.off()



table(tf_fam$Family)




# Write Tables for Fluidigm -----------------------------------------------

# ap2 top 1000 fluidigm

pcx %>% 
  arrange(PC1) %>%
  .[1:1000, ] %>%
  filter(Family == "AP2-EREBP") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "ap2 expressed in Branch Meristem") %>%
  write.csv2(., file = "../selected_genes/ap2-top-1000-pc1.csv")

# ap2 last 1000 fluidigm

pcx %>% 
  arrange(desc(PC1)) %>%
  .[1:1000, ] %>%
  filter(Family == "AP2-EREBP") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "ap2 expressed in Spikelet Meristem") %>%
  write.csv2(., file = "../selected_genes/ap2-last-1000-pc1.csv")


# MADS top 1000 fluidigm
pcx %>% 
  arrange(PC1) %>%
  .[1:1000, ] %>%
  filter(Family == "MADS") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "MADS expressed in Branch Meristem") %>%
  write.csv2(., file = "../selected_genes/MADS-top-1000-pc1.csv")

# zf-HD one gene - in branching meristem of african

pcx %>% 
  arrange(PC1) %>%
  # .[1:1000, ] %>%
  filter(locus_id == "LOC_Os11g03420") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "Expressed in Branch meristem of African species") %>%
  write.csv2(., file = "../selected_genes/zf−HD-top-1000-pc1.csv")

# DUF640 Expressed in branching, uncharacterized, shall we consider it?

pcx %>% 
  arrange(PC1) %>%
  .[1:1000, ] %>%
  left_join(mapman) %>%
  filter(grepl("DUF640", DESCRIPTION)) %>%
  select(locus_id, DESCRIPTION) %>%
  mutate(why = "Expressed in branching, uncharacterized, shall we consider it?") %>%
  write.csv2(., file = "../selected_genes/DUF640-top-1000-pc1.csv")
 
# LOC_Os03g08500LOC_Os04g23550 −−OsNAC19 down in glaberrima

# how many NAC?
# NAC 19 changes root architecture
# also OsNAC6, root architecture

pcx %>% 
  arrange(PC1) %>%
  .[1:1000, ] %>%
  filter(Family == "NAC") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "NAC expressed in Branch Meristem - they influence root architecture?") %>%
  write.csv2(., file = "../selected_genes/NAC-top-1000-pc1.csv")


### HOMEOBOX - but the list might be uncomplete
# HB top 1000 fluidigm

pcx %>% 
  arrange(PC1) %>%
  .[1:1000, ] %>%
  filter(Family == "HB") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "HB expressed in Branch Meristem") %>%
  write.csv2(., file = "../selected_genes/HB-top-1000-pc1-MAYBE-SOME-MISSING.csv")

# ap2 last 1000 fluidigm

pcx %>% 
  arrange(desc(PC1)) %>%
  .[1:1000, ] %>%
  filter(Family == "HB") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "HB expressed in Spikelet Meristem") %>%
  write.csv2(., file = "../selected_genes/HB-last-1000-pc1-MAYBE-SOME-MISSING.csv")

# MADS top 1000 fluidigm PC3
pcx %>% 
  arrange(PC3) %>%
  .[1:1000, ] %>%
  filter(Family == "MADS") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "MADS strange behaviour") %>%
  write.csv2(., file = "../selected_genes/MADS-top-1000-pc3.csv")

# and last
pcx %>% 
  arrange(desc(PC3)) %>%
  .[1:1000, ] %>%
  filter(Family == "MADS") %>%
  left_join(annos) %>%
  select(locus_id, symbol) %>%
  mutate(why = "MADS strange behaviour") %>%
  write.csv2(., file = "../selected_genes/MADS-last-1000-pc3.csv")
