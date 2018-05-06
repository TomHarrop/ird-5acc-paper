
# Set up environment ------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(DESeq2)
load("../data/rlog-pca.Rdata")
# Set an enstablished colour palette
color_palette <- c("blue", "goldenrod")

# One-dimentional plot PC5 (splits species) -------------------------------

fig_pc5 <- ggplot(pcx, aes(x = PC5, y = 0, 
                         colour = stage, 
                         label = ID)) +
  geom_hline(aes(yintercept=0), colour = "grey") +
  geom_point(size = 2, alpha = .5) +
  ggrepel::geom_text_repel() +
  scale_color_manual(values = color_palette) +
  theme_bw() +
  theme(line = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(expression(~5^th~Component))

svg("../fig/fig-pc5.svg",
    width = 7,
    height = 1.8)
fig_pc5
dev.off()


# Explore loadings of pc5 -------------------------------------------------


tst <- pcx %>% group_by(accession, stage) %>% nest(PC5) %>%
  mutate(range = map(.x = .$data, .f = range))
tst
tst$range

# 2-dimensional plots -----------------------------------------------------

# Plot PC

plot_pc <- function(pc_x, pc_y)
{
  ggplot(pcx, aes_string(x = pc_x, y = pc_y,
                         colour = "accession",
                         pch = "stage")) +
    geom_point(size = 3, alpha = .7) +
    theme_bw()
}

pdf("../fig/fig-pc-all.pdf",
    width = 7,
    height = 6)
# Define a set of alternate pairs of components that must be plotted
to_plot <- tibble(pc_x = paste0("PC", (1:7)*2 - 1),
                  pc_y = paste0("PC", (1:7)*2))
# And plot them all
to_plot %>% pmap(plot_pc)
dev.off()

# One-dimentional plot (PC5 splits species) -------------------------------

pcx <- pcx %>%
  gather(PC1:PC30, key = pc_id, value = pcx) %>%
  filter(pc_id %in% paste0("PC", 1:8))

fig_pc <- ggplot(pcx, aes(x = pcx, y = 0,
                          colour = accession,
                          pch = stage,
                          label = ID)) +
  geom_hline(aes(yintercept=0), colour = "grey") +
  geom_point(size = 2, alpha = .5) +
  ggrepel::geom_text_repel() +
  # scale_color_manual(values = color_palette) +
  facet_wrap(facets = "pc_id",
             scales = "free_x",
             ncol = 1) +
  theme_bw() +
  theme(line = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  # ggtitle(expression(~5^th~Component))
  ggtitle("Main Components")


# svg("../fig/fig-pc5.svg",
#     width = 7,
#     height = 9)
# fig_pc
# dev.off()
