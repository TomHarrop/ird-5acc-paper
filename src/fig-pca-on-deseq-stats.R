library(tidyverse)
# library(DESeq2)

load("../data/all-sep-deseq.Rdata")

spc_res <- spc_res %>%
  map(~.[, c("locus_id", "stat")]) %>%
  map(as.data.frame) 

change_name <- function(i) {
  colnames(spc_res[[i]])[2] <- paste("stat", i, sep = "_")
  return(spc_res[[i]])
}

spc_res <- names(spc_res) %>% map(change_name) 

spc_res <- spc_res %>% reduce(full_join, by = "locus_id")

rownames(spc_res) <- spc_res$locus_id
spc_res$locus_id <- NULL
spc_res[is.na(spc_res)] <- 0

pc <- prcomp(spc_res, scale. = T)

pc$rotation

pc$sdev

library(fgsea)
source("helper-functions.R")
load("../data/mapman.Rdata")

pcx <- as.data.frame(pc$x)
pcx$locus_id <- rownames(pcx)

mapman_list <- c(split(mapman$IDENTIFIER, mapman$level1),
                 split(mapman$IDENTIFIER, mapman$level2),
                 split(mapman$IDENTIFIER, mapman$level3),
                 split(mapman$IDENTIFIER, mapman$level4),
                 split(mapman$IDENTIFIER, mapman$level4),
                 split(mapman$IDENTIFIER, mapman$level6),
                 split(mapman$IDENTIFIER, mapman$level7))

pdf("../fig/fig-pca-on-deseq-stats.pdf",
    height = 12, width = 20)
gsea_stats <- map(pcx[, 1:5],
                  ~test_gsea(set_names(., 
                                       nm = pcx$locus_id),
                             mapman = mapman_list,
                             plot_top = T))
dev.off()

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

cutoff <- 200

pc3 <- pcx  %>%
  select(PC3, locus_id)

plot_all <- function(locus_ids,
                     tag = "PC3",
                     height = cutoff/1.3) {
  p <- locus_ids %>%
    get_expression(dds) %>%
    left_join(annos) %>%
    left_join(tf_fam) %>%
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
  
  # pdf(file = paste0("../fig/fig-tmp-",
  pdf(file = paste0("../fig/fig-tmp-pc-stat-deseq-",
                    tag,
                    ".pdf"),
      height = height,
      width = 18)
  print(p)
  dev.off()
}

plot_all(pcx %>% arrange(desc(PC2)) %>% .$locus_id %>% .[1:cutoff])
