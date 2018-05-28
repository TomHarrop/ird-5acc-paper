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

pheno_cali %>%
  ggplot(., aes(x = PbIntL,
                y = pbn)) +
  geom_point(alpha = .2) +
  facet_wrap(facets = "Origin") +
  theme_bw()

pheno_cali %>%
  ggplot(., aes(x = RL,
                y = pbn)) +
  geom_point(alpha = .2) +
  # geom_hex() +
  facet_wrap(facets = "Origin") +
  theme_bw()

pheno_cali %>%
  ggplot(., aes(x = pbn, y = sbn)) +
  geom_point(alpha = .2) +
  facet_wrap(facets = "Origin")
