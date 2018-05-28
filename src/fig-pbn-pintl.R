library(tidyverse)
library(gridExtra)

load("../data/phenotypes.Rdata")

pheno_mnp %>%
  ggplot(., aes(x = species, y = pbn)) +
  # geom_jitter(width = .1, height = 0) + 
  geom_count() +
  theme_bw()


pheno_mnp %>%
  ggplot(., aes(x = `Rachis_length (RL)`,
                y = pbn)) +
  geom_point() +
  theme_bw()

pheno_mnp %>%
  ggplot(., aes(x = `primary branch internode average length (PBintL in cm)`,
                y = pbn)) +
  geom_point() +
  # facet_wrap(facets = "species") +
  theme_bw()

svg("../fig/fig-pbn-pintl.svg")
pheno_cali %>%
  ggplot(., aes(x = pbn,
                y = PbIntL)) +
  geom_point(alpha = .2) +
  # geom_boxplot() +
  facet_wrap(facets = "Origin") +
  ylab("Primary Branch Internode Length")+
  xlab("Primary Branch Number") +
  theme_bw()
dev.off()

svg("../fig/fig-RL-pbn.svg")
pheno_cali %>%
  ggplot(., aes(x = pbn,
                y = RL)) +
  geom_point(alpha = .2) +
  # geom_hex() +
  facet_wrap(facets = "Origin") +
  ylab("Rachis Length")+
  xlab("Primary Branch Number") +
  theme_bw()
dev.off()

svg("../fig/fig-pbn-sbn.svg")
pheno_cali %>%
  ggplot(., aes(x = pbn, y = sbn)) +
  geom_point(alpha = .2) +
  ylab("Secondary Branch Number") +
  xlab("Primary Branch Number") +
  facet_wrap(facets = "Origin") +
  theme_bw()
dev.off()

pheno_cali %>%
  ggplot(., aes(x = sbn, y = pbl)) +
  geom_point(alpha = .2) +
  # ylab("Secondary Branch Number") +
  # xlab("Primary Branch Number") +
  facet_wrap(facets = "Origin") +
  theme_bw()

pheno_cali$PbIntL %>%
  ggplot(., aes(x = sbn, y = )) +
  geom_point(alpha = .2) +
  # ylab("Secondary Branch Number") +
  # xlab("Primary Branch Number") +
  facet_wrap(facets = "Origin") +
  theme_bw()
