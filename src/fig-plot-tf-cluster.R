library(tidyverse)
library(DESeq2)

source("helper-functions.R")
color_palette <- c("blue", "goldenrod")

dds <- readRDS("../data-raw/dds.Rds")
load("../data/rlog-pca.Rdata")
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

plot_cluster <- function(cl) {
  p <- get_expression(locus_ids = cls[[cl]]$locus_id,
                      dds = dds) %>%
    left_join(annos) %>%
    left_join(tf_fam) %>%
    left_join(mapman) %>%
    mutate(locus_id = as_factor(locus_id)) %>%
    plot_norm_expr() +
    ggtitle(paste0("Cluster", cl)) +
    facet_wrap(facets = c("locus_id",
                          "symbol",
                          "Family",
                          "DESCRIPTION"),
               scales = "free_y",
               ncol = 5,
               labeller = label_wrap_gen(width = 50,
                                         multi_line = T))
  
  nplots <- length(levels(p$data$locus_id))
  pdf(file = paste0("../fig/fig-tmp-cluster-",
                    # svg(file = paste0("../fig/fig-tmp-",
                    cl,
                    "-all.pdf"),
      # ".svg"),
      height = ceiling(nplots/5)*4,
      width = ifelse(nplots < 5, nplots * 4, 20))
  print(p)
  dev.off()
}

cls <- read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  dplyr::rename(locus_id = "MsuID") %>%
  split(.$cluster)

names(cls) %>% walk(plot_cluster)


