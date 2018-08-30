library(tidyverse)
library(pheatmap)
library(ggpubr)
library(DESeq2)


source("../src/helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
load("../data/msu-to-rapdb.Rdata"); rap2msu <- dict; rm(dict)

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

pc5 <- pcro %>%
  select(PC5, locus_id)

# check also clusters
clusters <- read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  mutate(locus_id = MsuID) %>%
  select(locus_id, cluster) %>%
  mutate(cluster = paste0("cluster", cluster))

# AP2s From Nakano 2006 --------------------------------------------------------

ap2s <- "http://www.plantphysiol.org/highwire/filestream/120654/field_highwire_adjunct_files/1/73783supptable_III_rv_1202.xls"

dest <- "../data-raw/ap2s-nakano.xls"

if(!file.exists(dest)) ap2s %>% download.file(destfile = dest)

ap2s_nakano <- readxl::read_excel(path = dest, range = "A2:G141") %>%
  mutate(locus_id = paste0("LOC_", `TIGR Locus identifier`))

# AP2s from Sharoni 2010 ---------------------------------------------------------
 
dest <- "../data-raw/ap2s-sharoni.zip"

if(!file.exists(dest)) {
  download.file("https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/pcp/52/2/10.1093_pcp_pcq196/5/pcq196_Supplementary_Data.zip?Expires=2147483647&Signature=FxY0bI9F1ooCm6XoS176tuLR5awBqvoTvfRPXk3xNzoI6sRFbWm7gEWUpqEtJzNdXrNu0AflCLTgHc4GG-ucfSpEjSuoVRIu~5nGqH~HHNwS~OaMrHpTu4NV5ZbT~XOxlDxb7PZXCUnXy44PIqQGXPErUpuqaz7A2OBdCEnyMClPQWrWvYtQwbCZJdePi-PjvD-wEribXVjUvdeX0puyN8EZy~puufSGiHxGxTRo387GgVyU05Ev3Uf6YrAQ-zJaUnFLB3~Km9mQh76uV91onEm8BLBLRggX0YqhBW5YktudYsiGaj50HuPvTi9jUnVOg5JFQpP4KJGQVdX5MmUvsw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA",
                destfile = dest)
}

ap2_sharoni_path <- "../data-raw/ap2s-sharoni/Supplementary_Table_S1.xls"

if(!file.exists(ap2_sharoni_path)) {
  unzip(dest,
        files = "Supplementary_Table_S1.xls",
        exdir = "../data-raw/ap2s-sharoni/")
}


ap2_sharoni <- bind_cols(read_excel(ap2_sharoni_path,
                                     range = "A4:A166",
                                     col_names = "locus_id"),
                          read_excel(ap2_sharoni_path,
                                     range = "N4:O166",
                                     col_names = c("subgroub", "subfamily"))) %>%
  mutate(locus_id = paste0("LOC_", locus_id),
         subfam = paste(subfamily, subgroub, sep = " - ")) %>%
  left_join(pc5) %>%
  # filter(!is.na(PC5)) %>%
  arrange(desc(PC5))

tst <- ap2_sharoni %>%
  left_join(annos)


# plot ap2s with subfam ---------------------------------------------------
  

p <- ap2_sharoni$locus_id %>%
  get_expression(dds = dds) %>%
  left_join(annos) %>%
  left_join(tf_fam) %>%
  left_join(mapman) %>%
  left_join(ap2_sharoni) %>%
  left_join(clusters) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id",
                        "symbol",
                        "Family",
                        "DESCRIPTION",
                        "subfam", 
                        "cluster"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 50,
                                       multi_line = T))

nplots <- length(levels(p$data$locus_id))
pdf(file = paste0("../fig/ap2-pc5-with-subfam.pdf"),
    # ".svg"),
    height = ceiling(nplots/5)*4,
    width = ifelse(nplots < 5, nplots * 4, 20))
print(p)
dev.off()


# Check also clusters -----------------------------------------------------

# read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
#   mutate(locus_id = MsuID) %>% 
#   inner_join(ap2_sharoni) %>% 
#   select(locus_id, names, cluster, subfam) %>% 
#   View()
#   