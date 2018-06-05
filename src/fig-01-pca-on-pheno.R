library(tidyverse)
library(ggfortify)
library(grid)
library(gridExtra)
library(ggpubr)  

load("../data/phenotypes.Rdata")


# PCA ---------------------------------------------------------------------

pc <- pheno_cali %>%
  select_at(vars(RL:TbN)) %>%
  select(-PanL) %>%
   as.data.frame() %>%
  prcomp(.,
         scale. = T,
         center = T)

pdf(file = "../fig/fig-01-pca-on-pheno.pdf")
autoplot(pc, data = pheno_cali, 
                     colour = "Type",
                     shape = "Species",
                     loadings = TRUE,
                     loadings.label = TRUE) +
  # geom_point(alpha = .3) +
  labs(caption = str_wrap("Fig. 1: PC Analysis of the scaled and centered 
                    panicle phenotype dataset
                    The first component explains 40% of the variance and splits wild and
                    cultivated species.
                    The first component splits traits that are intuitively related to high
                    yield, such as branch number and branch length, from traits that 
                    intuitively should correlate inversely with yield, such as
                    internode lenght.",
                          width = 70)) +
  theme_bw() +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1))
  dev.off()
plot(pheno_cali$spn ~ pc$x[, 1])


# PCA old ----------------------------------------------------------------


# ggplot(as.data.frame(pc$x)  %>%
#          cbind(Origin = pheno_cali$Origin),
#        aes(x = PC1, y = PC2, colour = Origin)) +
#   geom_point() +
#   theme_bw()
# 
# ggplot(pc$rotation %>%
#          as.data.frame(.) %>%
#          rownames_to_column(),
#        aes(x = PC1, y = PC2, label = rowname)) +
#   geom_text() +
#   theme_bw()

# Correlation and linear model --------------------------------------------

summary(lm(pheno_cali$spn ~ pc$x[, 1]))
ggplot(pheno_cali, aes(x = sbn, y = spn)) +
  geom_hex(bins = 30) +
  facet_wrap(facets =  "Species") +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 
ggplot(pheno_cali, aes(x = pbn, y = spn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = PanL, y = spn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = PanL, y = pbn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species")
ggplot(pheno_cali, aes(x = pbn, y = sbn)) +
  geom_hex(bins = 20) +
  facet_wrap(facets =  "Species") +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 

pheno_cali %>%
  split(.$Species) %>%
  map(~cor(.$pbn, .$sbn))


pheno_cali %>%
  split(.$Species) %>%
  map(~c(range(.$pbn), range(.$sbn)))

cor(pheno_cali$spn, pheno_cali$sbn)

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

pheno_cali %>%
  mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
  ggplot(aes(x = Species, y = naxil)) + 
  geom_boxplot()


pheno_cali %>%
  mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
  ggplot(aes(x = spn, y = naxil, colour = Species)) + 
  geom_point(alpha = .5) +
  # facet_wrap(facets = "Species") +
  theme_bw()

pheno_cali %>%
  mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
  split(.$Species) %>%
  map(~mean(.$naxil))


pheno_mnp %>%
  mutate(naxil = (spn - pbn)/(pbn + sbn + `Tb_nb (TbN)`)) %>%
  ggplot(aes(x = spn, y = naxil, colour = Species)) + 
  geom_point(alpha = .5) +
  # facet_wrap(facets = "Species") +
  theme_bw()

pheno_mnp %>%
  mutate(naxil = (spn - pbn)/(pbn + sbn + `Tb_nb (TbN)`)) %>%
  split(.$Species) %>%
  map(~mean(.$naxil))
