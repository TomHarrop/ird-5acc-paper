library(tidyverse)
load("../data/all-sep-deseq.Rdata")
source("helper-functions.R")
# load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
library(fgsea)


mapman <- mapman %>%
  dplyr::rename(locus_id = "IDENTIFIER")


mapman_list <- c(
  split(mapman$locus_id, mapman$level1),
  split(mapman$locus_id, mapman$level2),
  split(mapman$locus_id, mapman$level3)
  # split(mapman$locus_id, mapman$level4),
  # split(mapman$locus_id, mapman$level5),
  # split(mapman$locus_id, mapman$level6),
  # split(mapman$locus_id, mapman$level7)
)

plot_gsea <- function(res) 
{
  ge <- test_gsea(rnk = set_names(res$stat,
                                  nm = res$locus_id), 
                  mapman_list = mapman_list,
                  plot_top = T)
  return(ge)
}

pdf(file = "../fig/fig-explore-all-sep-deseq.pdf",
    width = 10, height = 8)
spc_res %>% map(plot_gsea)
dev.off()
