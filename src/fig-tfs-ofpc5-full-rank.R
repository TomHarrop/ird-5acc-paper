load("../data/rlog-pca.Rdata")
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

pcro <- pcro %>% 
  # as.data.frame(.) %>%
  # rownames_to_column() %>%
  # dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam) %>%
  left_join(mapman) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))

# Plot selected -----------------------------------------------------------

plot_kword <- function(kword = "AP2-EREBP") {
  p <- pcro %>%
    arrange(PC5) %>%
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

plot_mapman <- function(kword = "zinc finger") {
  p <- pcro %>%
    arrange(PC5) %>%
    filter(str_detect(DESCRIPTION, kword)) %>%
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

pdf("../fig/fig_ap2_from_pc5.pdf",
    width = 20,
    height = 60)
print(plot_kword())
dev.off()

pdf("../fig/fig_mads_from_pc5.pdf",
    width = 20,
    height = 40)
print(plot_kword("MADS"))
dev.off()

pdf("../fig/fig_zf_hd_from_pc5.pdf",
    width = 20,
    height = 25)
print(plot_kword("zf-HD"))
dev.off()


pdf("../fig/fig_NAC_from_pc5.pdf",
    width = 20,
    height = 30)
print(plot_kword("NAC"))
dev.off()

pdf("../fig/fig_HB_from_pc5.pdf",
    width = 20,
    height = 60)
print(plot_kword("HB"))
dev.off()

pdf("../fig/fig_SPB_from_pc5.pdf",
    width = 20,
    height = 20)
print(plot_kword("SBP"))
dev.off()

pdf("../fig/fig_zinc_finger_from_pc5.pdf",
    width = 20,
    height = 500)
print(plot_mapman())
dev.off()



table(tf_fam$Family)
