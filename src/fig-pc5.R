
# Set up environment ------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(DESeq2)
load("../data/pca-rlog.Rdata")
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
  theme(line = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(expression(~5^th~Component))

svg("../fig/fig-pc5.svg",
    width = 7,
    height = 1.8)
fig_pc5
dev.off()
