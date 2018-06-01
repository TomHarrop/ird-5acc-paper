library(tidyverse)
library(ggfortify)
library(gridExtra)

load("../data/phenotypes.Rdata")

pc <- pheno_cali %>%
  select_at(vars(RL:TbN)) %>%
  as.data.frame() %>%
  prcomp(.,
         scale. = T,
         center = T)

tst <- pheno_cali %>% 
  select(Bar_Code, Species) %>%
  distinct() %>%
  split(.$Species)

table(tst$Bar_Code, tst$Species)

ggplot(as.data.frame(pc$x)  %>%
         cbind(Origin = pheno_cali$Origin),
       aes(x = PC1, y = PC2, colour = Origin)) +
  geom_point() +
  theme_bw()

ggplot(pc$rotation %>%
         as.data.frame(.) %>%
         rownames_to_column(),
       aes(x = PC1, y = PC2, label = rowname)) +
  geom_text() +
  theme_bw()

summary(lm(pheno_cali$spn ~ pc$x[, 1]))
plot(pheno_cali$spn ~ pc$x[, 1])
ggplot(pheno_cali, aes(x = pbn, y = spn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = PanL, y = spn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = PanL, y = pbn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = sbn, y = spn)) +
  geom_hex(bins = 30) +
  facet_wrap(facets =  "Species") 
ggplot(pheno_cali, aes(x = pbn, y = sbn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")


plot(pheno_cali$spn ~ pheno_cali$pbn)
plot(pheno_cali$spn ~ pheno_cali$sbn)

pheno_mnp <- pheno_mnp %>%
  filter(complete.cases(.))

pc_mnp <- pheno_mnp %>%
  select_at(vars(`Rachis_length (RL)`:`secondary branch internode average length (SbintL in cm)`)) %>%
  as.data.frame() %>%
  # filter(complete.cases(.)) %>%
  prcomp(.,
         scale. = T,
         center = T)

ggplot(as.data.frame(pc_mnp$x)  %>%
         cbind(Origin = pheno_mnp$Origin),
       aes(x = PC1, y = PC2, colour = Origin)) +
  geom_point() +
  theme_bw()

pc_mnp$sdev

plot(pheno_mnp$spn ~ pc_mnp$x[,1])
plot(pheno_mnp$spn ~ pheno_mnp$sbn)
plot(pheno_mnp$spn ~ pheno_mnp$pbn)

