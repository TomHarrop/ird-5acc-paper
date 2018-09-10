library(tidyverse)
library(readxl)

dat <- read_excel(path = "../data-raw/Mutant_AP2_PhenotypingData.xlsx")

dat <- dat %>% 
  mutate_at(vars(RL:SpN), as.numeric) %>% 
  gather(RL:SpN, key = "measure", value = "value")


# Plot every mutant -------------------------------------------------------

# pdf("figures/mutants-boxplot.pdf", height = 14)
ggplot(dat, aes(x = `Target Gene`, y = value)) +
  geom_boxplot() +
  facet_grid(measure ~ Accession, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

