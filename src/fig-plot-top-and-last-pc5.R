library(tidyverse)

# Load datasets -----------------------------------------------------------

source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
load("../data/msu2rap.Rdata")


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


# Get pc5 -----------------------------------------------------------------

pc5 <- pcro %>%
  select(PC5, locus_id)

last_pc5 <- pc5 %>%
  top_n(-100, wt = PC5) %>%
  pull(locus_id) 

top_pc5 <- pc5 %>%
  top_n(100, wt = PC5) %>%
  pull(locus_id)


# Plot --------------------------------------------------------------------
plot_ids <- function(ids, file_path) 
{
  p <- ids %>%
    get_expression(dds = dds) %>%
    left_join(annos) %>%
    left_join(tf_fam) %>%
    left_join(mapman) %>%
    left_join(msu2rap) %>%
    mutate(locus_id = as_factor(locus_id)) %>%
    plot_norm_expr() +
    facet_wrap(facets = c("locus_id",
                          "symbol",
                          "Family",
                          "DESCRIPTION",
                          "rap_id"),
               scales = "free_y",
               ncol = 5,
               labeller = label_wrap_gen(width = 50,
                                         multi_line = T))
  
  nplots <- length(levels(p$data$locus_id))
  
  pdf(file = file_path,
      height = ceiling(nplots/5)*7,
      width = ifelse(nplots < 5, nplots * 4, 20))
  print(p)
  dev.off()
}

plot_ids(last_pc5, "../fig/fig-last-pc5.pdf")
plot_ids(top_pc5, "../fig/fig-top-pc5.pdf")
