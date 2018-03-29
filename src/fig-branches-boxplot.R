library(tidyverse)
library(gridExtra)

load("../data/phenotype-tidy.Rdata")


# Plot pbn sbn and spn together -------------------------------------------

xlab <- "species"

pbn <- ggplot(pheno_mnp, aes(x = species_simple, y = pbn)) +
  geom_boxplot() +
  geom_count(alpha = .7) +
  ylab("primary branch number") +
  xlab(xlab) +
  ggtitle("Primary Branches") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

sbn <- ggplot(pheno_mnp, aes(x = species_simple, y = sbn)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = .6, height = 0, width = .2) +
  ylab("secondary branch number") +
  xlab(xlab) +
  ggtitle("Secondary Branches") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

spn <- ggplot(pheno_mnp, aes(x = species_simple, y = spn)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = .6, height = 0, width = .2) +
  ylab("spikelet number") +
  xlab(xlab) +
  ggtitle("Spikelets") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

svg(file = "../fig/fig-branches-boxplot.svg",
    width = 7, height = 4)
grid.arrange(grobs = list(pbn, sbn, spn),
             layout_matrix = matrix(c(1,2,3), nrow = 1))
dev.off()
