library(tidyverse)
load("../data/rlog-pca.Rdata")
color_palette <- c("blue", "goldenrod")

pcx <- pcx %>%
  select_at(vars(accession:PC5)) %>%
  gather(PC1:PC5, key = "PC_id", value = "pc_x")

pdf("../fig/fig-04-PCA-rlog-VISUAL.pdf")
ggplot(pcx,
       aes(x = ID,
           fill = stage,
           # colour = stage,
           y = pc_x)) +
  geom_col() +
  # geom_point(size = 4) +
  geom_hline(yintercept = 0) +
  facet_grid(PC_id ~ accession, scales = "free") +
  theme_bw() +
  # scale_color_manual(values = color_palette) 
  scale_fill_manual(values = color_palette)
dev.off()
