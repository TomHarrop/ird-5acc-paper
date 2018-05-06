library(tidyverse)
library(pheatmap)
library(ggpubr)
color_palette <- c("blue", "goldenrod")

source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")
tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

pc7 <- pcro %>%
  select(PC7, locus_id)

prepare_for_plot <- function(dats) {
  dats <- get_expression(locus_ids = dats$locus_id,
                         dds = dds)
  
  dats <- dats %>%
    group_by(locus_id, species, stage) %>%
    # mutate(scaled_expression = scale(`normalized expression`)) %>%
    summarise(by_locus_species = median(by_locus_species)) %>%
    ungroup() %>%
    mutate(stage_species = paste(stage, species, sep = "_")) %>%
    select(locus_id, stage_species, by_locus_species) %>%
    spread(key = stage_species, value = by_locus_species) %>%
    left_join(annos) %>%
    left_join(tf_fam) %>%
    mutate(locus_name = paste(locus_id,
                              symbol,
                              Family,
                              # oryzr::LocToGeneName(.$locus_id)$symbols,
                              sep = "-")) %>%
    # select(-symbol) %>%
    distinct() %>%
    as.data.frame(.)
  
  rownames(dats) <- dats$locus_name
  return(dats)
}


pc7 <- pc7 %>% arrange(PC7)
bot_pc7 <- pc7 %>%
  # like this they are ordered by ranking
  .[1:200, ]  %>%
  left_join(prepare_for_plot(pc7[1:200, ])); rownames(bot_pc7) <- bot_pc7$locus_name

pc7 <- pc7 %>% arrange(desc(PC7))
top_pc7 <- pc7 %>% 
  # like this they are ordered by ranking
  .[1:200, ] %>%
  left_join(prepare_for_plot(pc7[1:200, ])); rownames(top_pc7) <- top_pc7$locus_name

tp <- pheatmap(top_pc7[, 3:12],
               main = "pc7, top 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               # cutree_rows = cut_tp,
               cluster_cols = F,
               cluster_rows = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9,
               filename = "../fig/heatmap-dump.pdf")

bt <- pheatmap(bot_pc7[, 3:12], 
               main = "pc7, last 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               # cutree_rows = cut_bt,
               cluster_cols = F,
               cluster_rows = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = 9,
               filename = "../fig/heatmap-dump.pdf")

pdf("../fig/fig-genes-of-pc7-explore.pdf",
    width = 14, height = 27)
ggarrange(plotlist = list(tp[[4]],
                          bt[[4]]))
dev.off()

plot_family <- function(dats, family, height = 4) {
  p <- dats %>%
    filter(Family == family) %>%
    .$locus_id %>%
    get_expression(dds) %>%
    left_join(annos) %>%
    plot_norm_expr() +
    facet_wrap(facets = c("locus_id",
                          "symbol"),
               scales = "free_y",
               ncol = 5)
  
  pdf(file = paste0("../fig/fig-tmp-",
                    family, "-",
                    substitute(dats),
                    ".pdf"),
      height = height,
      width = 12)
  print(p)
  dev.off()
}

plot_family(bot_pc7, family = "AP2-EREBP")

plot_tops <- function(dats, height = 4) {
  p <- dats %>%
    .[1:20, ] %>%
    .$locus_id %>%
    get_expression(dds) %>%
    left_join(annos) %>%
    plot_norm_expr() +
    facet_wrap(facets = c("locus_id",
                          "symbol"),
               scales = "free_y",
               ncol = 5)
  
  pdf(file = paste0("../fig/fig-tmp-",
                    "top20-",
                    substitute(dats),
                    ".pdf"),
      height = height,
      width = 12)
  print(p)
  dev.off()
}

plot_tops(bot_pc7, height = 12)
