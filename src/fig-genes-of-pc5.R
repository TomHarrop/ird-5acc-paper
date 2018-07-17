library(tidyverse)
library(pheatmap)
library(ggpubr)
color_palette <- c("blue", "goldenrod")

source("helper-functions.R")
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
    left_join(mapman) %>%
    mutate(locus_name = paste(locus_id,
                              symbol,
                              Family,
                              DESCRIPTION,
                              # oryzr::LocToGeneName(.$locus_id)$symbols,
                              sep = "-")) %>%
    # select(-symbol) %>%
    distinct() %>%
    as.data.frame(.)
  
  rownames(dats) <- dats$locus_name
  return(dats)
}

cutoff <- 400

pc5 <- pc5 %>% arrange(PC5)
bot_pc5 <- pc5 %>%
  # like this they are ordered by ranking
  .[1:cutoff, ]  %>%
  left_join(prepare_for_plot(pc5[1:cutoff, ])); rownames(bot_pc5) <- bot_pc5$locus_name

pc5 <- pc5 %>% arrange(desc(PC5))
top_pc5 <- pc5 %>% 
  # like this they are ordered by ranking
  .[1:cutoff, ] %>%
  left_join(prepare_for_plot(pc5[1:cutoff, ])); rownames(top_pc5) <- top_pc5$locus_name

tp <- pheatmap(top_pc5[, 3:12],
               main = "pc5, top 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               show_rownames = F,
               cutree_cols = 2,
               # cluster_cols = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = .2,
               filename = "../fig/heatmap-dump.pdf")

bt <- pheatmap(bot_pc5[, 3:12], 
               main = "pc5, last 200 genes",
               color = colorRampPalette(c( "white", "blue4"))(50),
               show_rownames = F,
               cutree_cols = 2,
               # cluster_cols = F,
               gaps_col = 5,
               cellwidth = 9,
               cellheight = .2,
               filename = "../fig/heatmap-dump.pdf")

svg("../fig/fig-genes-of-pc5.svg",
    width = 5.4, height = 3)
ggarrange(plotlist = list(tp[[4]],
                          bt[[4]]))
dev.off()


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

pdf("../fig/fig-genes-of-pc5-explore.pdf",
    # width = 14, height = 27)
    width = 27, height = 27)
ggarrange(plotlist = list(tp[[4]],
                          bt[[4]]))
dev.off()

# top_pc5_clust <- cbind(top_pc5,
#                        cluster = cutree(tp$tree_row, cut_tp))
# bot_pc5_clust <- cbind(bot_pc5, 
#                        cluster = cutree(bt$tree_row, cut_bt))
# 
# tmp <- get_expression(bot_pc5_clust %>%
#                         filter(cluster == 1) %>%
#                         .$locus_id,
#                       dds = dds) %>%
#   left_join(annos) %>%
#   mutate(species = factor(species,
#                           levels = c("japonica",
#                                      "barthii",
#                                      "glaberrima",
#                                      "rufipogon",
#                                      "indica")))
#   
# 
# pdf("../fig/tmp_genes_pc5.pdf",
#     height = 20,
#     width = 12)
# ggplot(tmp, aes(x = species,
#                 y = `normalized expression`,
#                 colour = stage,
#                 pch = domestication)) +
#   geom_point(size = 3, 
#              position = position_dodge(width = .5),
#              alpha = .8) +
#   scale_color_manual(values = color_palette) +
#   facet_wrap(facets = c("locus_id",
#                         "symbol"),
#              scales = "free_y",
#              ncol = 5) +
#   expand_limits(y=0) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   annotate("rect",
#            xmin = 1.5, xmax = 3.5,
#            ymin = -Inf, ymax = Inf,
#            alpha = .1)
# dev.off()
# 
# heat_by_clust <- function(dat, selection = "top") {
#   # tmp <- top_pc5_clust %>%
#   #   filter(cluster == 3)
#   if(nrow(dat) > 1) {
#     pheatmap(dat %>%
#                select_at(vars(PBM_barthii:SM_rufipogon)), #%>% t(.),
#              main = paste("pc5,", selection,  "200 genes - cluster", unique(dat$cluster)),
#              color = colorRampPalette(c( "white", "blue4"))(50),
#              labels_row = dat %>% .$locus_name,
#              # labels_col = oryzr::LocToGeneName(tmp$locus_id)$symbols,
#              cluster_cols = F,
#              gaps_col = 5,
#              cellwidth = 9,
#              cellheight = 9,
#              filename = paste0("../fig/fig-TMP-pc5-", selection, "200-cl",
#                               unique(dat$cluster),
#                               ".jpeg"))
#   }
# }
# top_pc5_clust %>% split(.$cluster) %>% walk(heat_by_clust)
# bot_pc5_clust %>% split(.$cluster) %>% walk(heat_by_clust, selection = "last")
# # tst <- oryzr::LocToGeneName(tmp$locus_id)$symbols

plot_family <- function(dats, family, height = 8) {
  p <- dats %>%
    # filter(Family == family) %>%
    filter(grepl(family, Family)) %>% #& !is.na(Family)) %>%
    .$locus_id %>%
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
  
  nplots <- length(levels(p$data$locus_id))
  pdf(file = paste0("../fig/fig-tmp-",
                    # svg(file = paste0("../fig/fig-tmp-",
                    family, "-",
                    substitute(dats),
                    ".pdf"),
                    # ".svg"),
      height = ceiling(nplots/5)*4,
      width = ifelse(nplots < 5, nplots * 4, 20))
  print(p)
  dev.off()
}


# top_pc5$locus_id %>%
#   pull_rap()

plot_all <- function(dats, height = cutoff/2) {
  
  
  p <- dats %>%
    # filter(Family == family) %>%
    # filter(grepl(family, Family)) %>% #& !is.na(Family)) %>%
    .$locus_id %>%
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
  pdf(file = paste0("../fig/fig-tmp-",
                    substitute(dats),
                    ".pdf"),
      height = height,
      width = 18)
  print(p)
  dev.off()
}


plot_mapman <- function(dats, mapman_word, height = 10) {
  p <- dats %>%
    # filter( == family) %>%
    filter(grepl(mapman_word, DESCRIPTION)) %>% #& !is.na(Family)) %>%
    .$locus_id %>%
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

  nplots <- length(levels(p$data$locus_id))
  # pdf(file = paste0("../fig/fig-tmp-",
  pdf(file = paste0("../fig/fig-tmp-",
                    mapman_word, "-",
                    substitute(dats),
                    ".pdf"),
      height = ceiling(nplots/5)*4,
      width = ifelse(nplots < 5, nplots * 4, 20))
  print(p)
  dev.off()
}


plot_family(bot_pc5, "AP2-EREBP")
plot_family(bot_pc5, "MADS")
plot_family(bot_pc5, "C2C2-YABBY")

plot_family(top_pc5, "AP2-EREBP")
plot_family(top_pc5, "MADS")
plot_family(top_pc5, "NAC")
plot_family(bot_pc5, "bHLH")
plot_family(top_pc5, "bHLH")
plot_family(top_pc5, "WRKY")
plot_family(top_pc5, "TRAF")


plot_all(top_pc5)
plot_all(bot_pc5)

plot_mapman(bot_pc5, "POEI", height = 6)
plot_mapman(bot_pc5, "LTPL", height = 6)
plot_mapman(bot_pc5, "ankyrin", height = 6)
plot_mapman(bot_pc5, "SPL", height = 6)
plot_mapman(bot_pc5, "zinc", height = 6)
plot_mapman(bot_pc5, "HD", height = 6)
plot_mapman(bot_pc5, "TCP", height = 6)
plot_mapman(bot_pc5, "YUC", height = 6)
plot_mapman(bot_pc5, "[Aa]uxin", height = 6)


plot_mapman(top_pc5, "zinc", height = 14)
plot_mapman(top_pc5, "TFL", height = 6)
plot_mapman(top_pc5, "LTPL", height = 6)
plot_mapman(top_pc5, "MYB", height = 6)
plot_mapman(top_pc5, "[Aa]uxin", height = 6)
plot_mapman(top_pc5, "[Gg]ibberellin", height = 6)
plot_mapman(top_pc5, "B3", height = 6)
plot_mapman(top_pc5, "DNA", height = 12)


