library(tidyverse)
library(ggfortify)
library(DESeq2)

load("../data/03-rld.Rdata")
colData(rld)

rlog_counts <- assay(rld)

pc <- prcomp(t(rlog_counts ), center = TRUE)

pcx <- as.data.frame(pc$x)
pcx$ID <- rownames(pcx)
pcx <- as_tibble(pcx)

samples <- colData(rld)
samples$ID <- rownames(samples)
samples <- as_tibble(as.data.frame(samples))

pcx <- inner_join(samples, pcx, by = "ID")

new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1),
                                         unit = "lines",
                                         valid.unit = 3L,
                                         class = "unit")

fig_pc5 <- ggplot(pcx, aes(x = PC5, y = 0, 
                         colour = stage, 
                         label = ID)) +
  geom_hline(aes(yintercept=0), colour = "grey") +
  geom_point(size = 2, alpha = .5) +
  ggrepel::geom_text_repel() +
  theme(line = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(expression(~5^th~Component))

svg("../fig/fig-pc5.svg",
    width = 7,
    height = 1.8)
fig_pc5

dev.off()


