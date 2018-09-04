library(tidyverse)
library(gridExtra)
load("../data/phenotypes.Rdata")

pheno_mnp <- pheno_mnp %>%
  mutate(Species = factor(Species,
                          levels = c("Oryza sativa japonica temp",
                                     "Ozyza Barthii",
                                     "Oryza glaberrima",
                                     "Oryza rufipogon",
                                     "Oryza sativa indica")))

# Plot --------------------------------------------------------------------

pbn <- ggplot(pheno_mnp, aes(x = Species, y = pbn)) +
  geom_boxplot() +
  geom_count(alpha = .7) +
  ylab("primary branch number") +
  ggtitle("Primary Branches") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

sbn <- ggplot(pheno_mnp, aes(x = Species, y = sbn)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = .6, height = 0, width = .2) +
  ylab("secondary branch number") +
  ggtitle("Secondary Branches") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

spn <- ggplot(pheno_mnp, aes(x = Species, y = spn)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = .6, height = 0, width = .2) +
  ylab("spikelet number") +
  ggtitle("Spikelets") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("rect",
           xmin = 1.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           alpha = .1)

pdf(file = "../fig/suppl-fig-01-phenotype-selected-acc-MNP.pdf",
    width = 7, height = 6)
grid.arrange(grobs = list(pbn, sbn, spn),
             layout_matrix = matrix(c(1,2,3), nrow = 1))
dev.off()
