library(tidyverse)
library(ggfortify)
library(grid)
library(gridExtra)
library(ggpubr)  

load("../data/phenotypes.Rdata")
color_palette <- c("#DAA520", "#0000F9",
                   "#FFE224", "#008AF9",
                   "#494949")

# PCA ---------------------------------------------------------------------

pc <- pheno_cali %>%
  select_at(vars(RL:spn)) %>%
  select(-PanL) %>%
   as.data.frame() %>%
  prcomp(.,
         scale. = T,
         center = T)
p <- autoplot(pc, data = pheno_cali, 
              loadings = TRUE,
              loadings.label = TRUE) 
# Get PC1 and  2 ready for plot -------------------------------------------


# pcx <- pc$x %>%
#   as.data.frame() %>%
#   bind_cols(pheno_cali) %>%
#   select(PC1, PC2, Species, Type, Species)

### the annotations from autoplot are in $layers[[3]]$data
# 
# pcx <- fortify(pc, data = pheno_cali) 
# 
# pcro <- pc$rotation %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   select(rowname, PC1, PC2) 
# %>%
#   mutate(PC1 = PC1*10, 
#          PC2 = PC2*7)

p_load <- p$layers[[3]]$data

p_1 <- p$data %>%
  ggplot(aes(x = PC1,
             y = PC2,
             colour = Species,
             pch = Species)) +
  geom_point(size = 2, alpha = .5) +
  annotate("segment",
           x = 0, xend = p_load$PC1,
           y = 0, yend = p_load$PC2,
           colour = color_palette[5],
           size = .4, alpha = 1,
           arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text",
           x = p_load$PC1*1.1,
           y = p_load$PC2*1.1, 
           label = p_load$rownames,
           colour = color_palette[5]) +
  labs(caption = str_wrap("Fig. 1: PC Analysis of the scaled and centered 
                    panicle phenotype dataset
                    The first component explains almost half of the variance and splits wild and
                    cultivated species.
                    The first component splits traits that are intuitively related to high
                    spikelet number, such as branch number and branch length, from traits that 
                    intuitively should correlate inversely with yield, such as
                    internode lenght. The highest loading given to secondary branch number and
                    spikelet number, which almost overlap.",
                          width = 70)) +
  theme_bw() +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1)) +
  scale_color_manual(values = color_palette)
pdf(file = "../fig/fig-01-pca-on-pheno.pdf")
p_1
dev.off()

autoplot(pc, data = pheno_cali, 
         colour = "Type",
         shape = "Species",
         loadings = TRUE,
         loadings.label = TRUE,
         loadings.colour = "black",
         loadings.label.colour = "black") +
  # geom_point(alpha = .3) +
  labs(caption = str_wrap("Fig. 1: PC Analysis of the scaled and centered 
                    panicle phenotype dataset
                    The first component explains almost half of the variance and splits wild and
                    cultivated species.
                    The first component splits traits that are intuitively related to high
                    spikelet number, such as branch number and branch length, from traits that 
                    intuitively should correlate inversely with yield, such as
                    internode lenght. The highest loading given to secondary branch number and
                    spikelet number, which almost overlap.",
                          width = 70)) +
  theme_bw() +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1)) +
  scale_color_manual(values = color_palette)
dev.off()

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
