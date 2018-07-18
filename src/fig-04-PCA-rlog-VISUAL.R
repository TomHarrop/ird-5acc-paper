library(tidyverse)
load("../data/rlog-pca.Rdata")
color_palette <- c("blue", "goldenrod")

# Only the first 5 PC are interesting? ------------------------------------

pcx <- pcx %>%
  select_at(vars(accession:PC5)) %>%
  gather(PC1:PC5, key = "PC_id", value = "pc_x")


# Plot and save -----------------------------------------------------------

p <- ggplot(pcx,
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


pdf("../fig/fig-04-PCA-rlog-VISUAL.pdf")
print(p)
dev.off()
