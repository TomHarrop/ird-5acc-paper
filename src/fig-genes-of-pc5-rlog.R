library(tidyverse)
library(pheatmap)
library(ggpubr)
library(DESeq2)

color_palette <- c("blue", "goldenrod")

source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

pc5 <- pcro %>%
  select(PC5, locus_id)

locus_ids <- "LOC_Os01g53370"
get_expr_rld <- function(locus_ids, rld, mapman = NULL) {
  
  # filter locus_ids
  excluded <- locus_ids[! locus_ids %in% rownames(rld)]
  if(length(excluded) > 0) { 
    warning(paste("These locus ids are not contained in the rld dataset:",
                  paste(excluded, collapse = ", ")))
  }
  
  locus_ids <- locus_ids[locus_ids %in% rownames(rld)]
  # first subset the DeSeqdataset, it makes everything
  # Downstream lighter
  rld <- rld[locus_ids, ]
  # print(rld)
  
  # We decided to variance stabilizing normalized counts
  dat <- as.data.frame(SummarizedExperiment::assay(rld))
  
  # Use locus id as rownames
  dat$locus_id <- rownames(dat)
  
  
  # get sample annotations
  cdata <- as.data.frame(SummarizedExperiment::colData(rld))
  cdata$sample <- rownames(cdata)
  
  # wrangle and merge
  dat <- dat %>%
    gather(colnames(rld), key = "sample", value = "normalized expression") %>%
    left_join(cdata, by = "sample") %>%
    dplyr::rename(species = accession)
  # str(dat)
  
  # merge mapman ids
  if(!is.null(mapman)) {
    mapman <- mapman %>% 
      dplyr::rename(locus_id = IDENTIFIER) %>%
      select(locus_id, DESCRIPTION) %>%
      distinct()
    # str(mapman)
    dat <- dat %>%
      left_join(mapman, by = "locus_id")
  }
  
  # scale by locus
  dat <- dat %>%
    group_by(locus_id) %>%
    mutate(by_locus = scale(`normalized expression`,
                            center = TRUE, scale = FALSE)) %>%
    ungroup()
  
  # scale by locus and species
  dat <- dat %>%
    group_by(species, locus_id) %>%
    mutate(by_locus_species = scale(`normalized expression`,
                                    center = TRUE, scale = FALSE)) %>%
    ungroup()
  
  return(dat)
}

prepare_for_plot <- function(dats) {
  dats <- get_expr_rld(locus_ids = dats$locus_id,
                         rld = rld)
  
  dats <- dats %>%
    group_by(locus_id, species, stage) %>%
    # mutate(scaled_expression = scale(`normalized expression`)) %>%
    summarise(by_locus_species = median(by_locus_species)) %>%
    ungroup() %>%
    mutate(stage_species = paste(stage, species, sep = "_")) %>%
    select(locus_id, stage_species, by_locus_species) %>%
    spread(key = stage_species, value = by_locus_species) %>%
    left_join(annos) %>%
    mutate(locus_name = paste(locus_id,
                              symbol,
                              oryzr::LocToGeneName(.$locus_id)$symbols,
                              sep = "-")) %>%
    # select(-symbol) %>%
    distinct() %>%
    as.data.frame(.)
  
  rownames(dats) <- dats$locus_name
  return(dats)
}

pc5 <- pc5 %>% arrange(PC5)
bot_pc5 <- pc5 %>% 
  .[1:200, ]  %>%
  left_join(prepare_for_plot(pc5[1:200, ]))

pc5 <- pc5 %>% arrange(desc(PC5))
top_pc5 <- pc5 %>% 
  .[1:200, ] %>%
  left_join(prepare_for_plot(pc5[1:200, ]))

cut_tp <- 9
cut_bt <- 6
tp <- pheatmap(top_pc5[, 3:12],
               main = "pc5, top 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               # cutree_rows = cut_tp,
               cluster_cols = F,
               cluster_rows = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9,
               filename = "../fig/heatmap-dump.pdf")

bt <- pheatmap(bot_pc5[, 3:12], 
               main = "pc5, last 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               # cutree_rows = cut_bt,
               cluster_cols = F,
               cluster_rows = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9,
               filename = "../fig/heatmap-dump.pdf")

pdf("../fig/fig-genes-of-pc5-rlog.pdf",
    width = 12, height = 27)
ggarrange(plotlist = list(tp[[4]],
                          bt[[4]]))
dev.off()

# top_pc5_clust <- cbind(top_pc5,
#                        cluster = cutree(tp$tree_row, cut_tp))
# bot_pc5_clust <- cbind(bot_pc5, 
#                        cluster = cutree(bt$tree_row, cut_bt))
