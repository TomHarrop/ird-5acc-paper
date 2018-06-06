# Plot result of z-score clustering


# Load Data ---------------------------------------------------------------


library(tidyverse)
library(readr)
dds <- readRDS("../data-raw/dds.Rds")
source("helper-functions.R")
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


clusts <- read_csv(file = "../data-raw/annotated_clusters_scaled_zstat.csv") %>% 
  split(.$cluster)


# Plot with annos ---------------------------------------------------------

plot_pdf_clust <- function(i) {
  pdf(paste0("../fig/fig-clust-stat-",
             i,
             ".pdf"),
      width = 18,
      height = nrow(clusts[[i]])/1.2)
  p <- get_expression(clusts[[i]]$MsuID,
                      dds = dds) %>%
    left_join(tf_fam)%>%
    left_join(annos) %>%
    left_join(mapman) %>%
    plot_norm_expr(.) +
    facet_wrap(facets = c("locus_id",
                          "symbol",
                          "Family",
                          "DESCRIPTION"),
               scales = "free_y",
               ncol = 5,
               labeller = label_wrap_gen(width = 50,
                                         multi_line = T))
  print(p)
  dev.off()
}

walk(names(clusts), plot_pdf_clust)

