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

# cl4 <- cls %>% filter(cluster == 4)
cl5 <- cls %>% filter(cluster == 5)

# Dot plot ---------------------------------------------------------------

dat <-  
  ap2_exp %>%
  filter(locus_id %in% cl5$locus_id) %>%
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj",
                                     "Ob", "Og",
                                     "Or", "Osi")))



# Plot -----------------------------------------------------


plts <- 
  dat %>% 
  {lineplot_fluidigm(nm = "",
                    dat = .,
                    alpha = .8)} +
  theme(strip.text = element_text(size = 8))

# +
#   theme(text = element_text(size = 8),
#         strip.text = element_text(size = 4))

# p <- cowplot::plot_grid(plts[[2]], plts[[3]],
#                         nrow = 2,
#                         rel_heights = c(3, 4) + 1.5,
#                         labels = c("A", "C")) %>%
#   cowplot::plot_grid(., plts[[1]],
#                      labels = c("", "B")) %>%
#   cowplot::add_sub(., str_wrap("qPCR confirms the behaviour of selected
#                                genes of cluster 4 and cluster 5.
#                                We have also measured gene expressio in two
#                                additional developmental stages. Stage 1 is a Rachis
#                                Meristem, Stage2 is a Branch Meristem, Stage 3
#                                is a Spikelet Meristem, Stage 4 is a Developing
#                                Spikelet."),
#                    size = 11) %>%
#   cowplot::ggdraw()


pdf("../fig/suppl-fig-fluidigm-ap2-hb.pdf",
    height = 9, width = 9,
    paper = "a4")
# do the inches really matter?
# PDF rescales (vector)
plts %>% print()
dev.off()


  