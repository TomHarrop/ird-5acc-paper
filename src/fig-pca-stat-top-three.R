library(tidyverse)
load("../data/all-sep-deseq.Rdata")

pcro <- pc_spc$rotation %>%
  as.data.frame(.) %>%
  mutate(species = rownames(.))

p <- ggplot(pcro %>%
              select_at(vars(PC1:PC3, species)) %>%
              gather(PC1:PC3, key = "Principal Component", value = "rotation"),
            aes(x = species, y = rotation)) + 
  geom_point(size = 3) +
  facet_grid(. ~ `Principal Component`) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

svg(filename = "../fig/fig-pca-stat-top-three.svg", 
    width = 7,
    height = 3)
print(p)
dev.off()
