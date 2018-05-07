library(tidyverse)
dds <- readRDS("../data-raw/dds.Rds")
load("../data/indica.Rdata")
source("helper-functions.R")
load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
library(fgsea)

mapman <- mapman %>%
  dplyr::rename(locus_id = "IDENTIFIER")

annos <- annos %>%
  dplyr::rename(locus_id = "MSU")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

tst <- tt_indica[tt_indica$padj < .05, ] %>%
# tst <- tt_indica[1:70, ] %>%
  left_join(annos) %>%
  left_join(tf_fam)


p <- get_expression(tst$locus_id, dds = dds) %>%
  left_join(annos) %>%
  left_join(mapman) %>% 
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id",
                        "symbol",
                        "DESCRIPTION"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 50))

pdf(file = "../fig/fig-explore-indica.pdf",
# svg(file = "../fig/fig-explore-indica.svg",
    height = 90,
    # height = 30,
    width = 17)
print(p)
dev.off()

# plot(tt_indica$stat, tt_indica$pvalue)

mapman_list <- c(
  split(mapman$locus_id, mapman$level1),
  split(mapman$locus_id, mapman$level2),
  split(mapman$locus_id, mapman$level3)
  # split(mapman$locus_id, mapman$level4),
  # split(mapman$locus_id, mapman$level5),
  # split(mapman$locus_id, mapman$level6),
  # split(mapman$locus_id, mapman$level7)
  )

ge <- test_gsea(rnk = set_names(tt_indica$stat,
                                nm = tt_indica$locus_id), 
                mapman_list = mapman_list,
                plot_top = T)
