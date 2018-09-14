library(tidyverse)
library(viridis)
library(pheatmap)
library(grid)
library(gridExtra)
source("../analyze-fluidigm/helper_functions_fluidigm.R")  


load("../data/ap2_fluidigm.R")
load("../data/hb_fluidigm.R")

cls <- read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  dplyr::rename(locus_id = "MsuID")

cl4 <- cls %>% filter(cluster == 4)
cl5 <- cls %>% filter(cluster == 5)

ap2_heat5 <- ap2_exp %>%
  filter(locus_id %in% cl5$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c5a <- pheatmap(mat = ap2_heat5 %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50), main = "AP2 cluster5",
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat5$locus_id)


ap2_heat4 <- ap2_exp %>%
  filter(locus_id %in% cl4$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c4a <- pheatmap(mat = ap2_heat4 %>% select(-locus_id),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50), main = "AP2 cluster4",
         cluster_cols = F,
         gaps_col = 1:4*5,
         cellwidth = 9,
         cellheight = 9,
         labels_row = ap2_heat4$locus_id)

hb_heat4 <- hb_exp %>%
  filter(locus_id %in% cl4$locus_id) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat()

c4h <- pheatmap(mat = hb_heat4 %>% select(-locus_id),
                # scale = "row",
                # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
                # color = colorRampPalette(c( "white", "blue4"))(50),
                color = viridis_pal()(50), main = "HB cluster4",
                cluster_cols = F,
                gaps_col = 1:4*5,
                cellwidth = 9,
                cellheight = 9,
                labels_row = hb_heat4$locus_id)

pdf("../fig/fig-06-fludigm-ap2-hb-DRAFT.pdf")
grid.arrange(c4h[[4]], c4a[[4]], c5a[[4]])
dev.off()


# Dot plot ---------------------------------------------------------------

dat <- list(hb_cl4 = hb_exp %>%
              filter(locus_id %in% cl4$locus_id),
            ap2_cl4 = ap2_exp %>%
              filter(locus_id %in% cl4$locus_id),
            ap2_cl5 = ap2_exp %>%
              filter(locus_id %in% cl5$locus_id)) %>%
  map(., ~scale_tidy_fluidigm(.) %>%
        mutate(species = factor(species,
                                levels = c("Osj",
                                           "Ob", "Og",
                                           "Or", "Osi"))))

dat <- map(names(dat), .f = ~mutate(dat[[.]],
                           new_var = .))

plts <- map(dat, ~ggplot(., aes(x = stage,
                                y = expression)) +
              geom_point(size = 2, col = "darkgrey") +
              geom_smooth(aes(x = stage %>%
                                as_factor(.) %>%
                                as.numeric(.)),
                          se = FALSE,
                          colour = "black") +
              facet_grid(target_name ~ species,
                         scales = "free_y") +
              # facet_grid(species ~ target_name, 
              #            scales = "free_x") +
              theme_bw() +
              theme(axis.text.x = element_text(hjust = 0,
                                               vjust = .5,
                                               angle = 270)) +
              labs(title = unique(.$new_var),
                   y = "Relative Expression")
            ) %>%
  map(., ggplotGrob)

pdf("../fig/fig-06-fluidigm-ap2-hb-dotline.pdf",
    height = 14)
    # width = 12)
g <- rbind(plts[[1]], plts[[2]], plts[[3]],
           size = "first")
# g <- cbind(plts[[1]], plts[[2]], plts[[3]],
#            size = "first")
grid.draw(g)
dev.off()


# Test coord fixed --------------------------------------------------------


plts_2 <- map(dat, ~ggplot(., aes(x = stage,
                                y = expression)) +
              geom_point(size = 2, col = "darkgrey") +
              geom_smooth(aes(x = stage %>%
                                as_factor(.) %>%
                                as.numeric(.)),
                          se = FALSE,
                          colour = "black") +
              facet_grid(target_name ~ species,
                         scales = "free_y", ) +
              # coord_fixed() +  
              # facet_grid(species ~ target_name, 
              #            scales = "free_x") +
              theme_bw() +
              theme(axis.text.x = element_text(hjust = 0,
                                               vjust = .5,
                                               angle = 270)) +
              labs(title = unique(.$new_var),
                   y = "Relative Expression")
) 

ph  <- 4

pdf("../fig/fig-05-fluidigm-ap2-hb-TEST-ratio.pdf",
    height = 8)
# width = 12)
grid.arrange(plts_2[[1]],
             plts_2[[2]],
             plts_2[[3]], 
             layout_matrix = rbind(c(rep(1, 4), rep(2, 4)),
                                   c(rep(1, 4), rep(2, 4)),
                                   c(rep(1, 4), rep(2, 4)),
                                   c(rep(1, 4), rep(3, 4)),
                                   c(rep(1, 4), rep(3, 4)),
                                   c(rep(1, 4), rep(3, 4)),
                                   c(rep(1, 4), rep(3, 4))))
# g <- cbind(plts[[1]], plts[[2]], plts[[3]],
#            size = "first")
# print(p)
dev.off()


# ggpubr::ggarrange(plts[[1]], plts[[2]], plts[[3]], ncol = 3)
# 
# # grid.arrange(plts[[1]], plts[[2]], plts[[3]])
# 
# grid.draw(gtable:::rbind_gtable(plts[[1]],
#                                 plts[[2]], "first"))
# 
# 
# pdf("../fig/fig-06-TEST.pdf", height = 20)
# grid.arrange(grobs = lapply(
#   list(plts[[1]],
#        plts[[2]],
#        plts[[3]]),
#   egg::set_panel_size,
#   width = unit(2, "cm"),
#   height = unit(1, "in")
# ))
# dev.off()
# 
# library(grid)
# library(gtable)
# g1 <- ggplotGrob(plts[[1]])
# g2 <- ggplotGrob(plts[[2]])
# g <- cbind(g1, g2, size = "first")
# 
# tst <- hb_exp %>%
#   filter(locus_id %in% cl4$locus_id) %>%
#   scale_tidy_fluidigm() %>%
#   mutate(species = factor(species,
#                           levels = c("Osj",
#                                      "Ob", "Og",
#                                      "Or", "Osi")))
# 
# 
#   ggplot(tst, aes(x = stage,
#              y = expression)) +
#   geom_point(size = 2, col = "darkgrey") +
#   geom_smooth(aes(x = stage %>%
#                     as_factor(.) %>%
#                     as.numeric(.)),
#               se = FALSE,
#               colour = "black") +
#   facet_grid(target_name ~ species, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(hjust = 0,
#                                    vjust = .5,
#                                    angle = 270))
# 

  