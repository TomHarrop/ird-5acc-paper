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


# Vertical arrangement ----------------------------------------------------


plts <- names(dat) %>%
  map(., ~lineplot_fluidigm(nm = .,
                            dat = dat[[.]],
                            alpha = .8))


pdf("../fig/fig-05-fluidigm-ap2-hb-dotline.pdf",
    height = 14)

g <- plts %>%
  map(., ggplotGrob) 
rbind(g[[1]], g[[2]], g[[3]],
           size = "first") %>%
  grid.draw()

dev.off()

ph  <- 4


# two columns -------------------------------------------------------------


pdf("../fig/fig-05-fluidigm-ap2-hb-TEST-ratio.pdf",
    height = 8)
cowplot::plot_grid(plts[[2]], plts[[3]],
                   nrow = 2,
                   rel_heights = c(3, 4) + 1.5,
                   labels = c("A", "C")) %>%
  cowplot::plot_grid(., plts[[1]],
                     labels = c("", "B")) %>%
  cowplot::add_sub(., str_wrap("")) %>%
  cowplot::ggdraw()
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

  