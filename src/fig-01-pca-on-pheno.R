library(tidyverse)
library(ggfortify)
library(grid)
library(gridExtra)
library(ggpubr)  
library(viridis)

load("../data/phenotypes.Rdata")
# color_palette <- viridis(20)[c(7, 19, 8, 20)]
color_palette <- viridis(20)[c(6, 14, 8, 16)]
# color_palette <- cividis(20)[c(1, 11, 2, 12)]
p_ch <- c(22, 24, 23, 25)


# Adjust species labels ---------------------------------------------------

pheno_cali <- pheno_cali %>%
  mutate(Species = paste0("O. ", Species))

# PCA ---------------------------------------------------------------------

pc <- pheno_cali %>%
  select_at(vars(RL:spn)) %>%
  select(-PanL) %>%
   as.data.frame() %>%
  prcomp(.,
         scale. = T,
         center = T)

var_expl <- summary(pc)$importance[2,]*100

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
             fill = Species,
             colour = Species,
             pch = Species)) +
  geom_point(size = 1,
             alpha = .5) +
  annotate("segment",
           x = 0, xend = p_load$PC1,
           y = 0, yend = p_load$PC2,
           colour = "#494949",
           size = .4, alpha = 1,
           arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text",
           x = p_load$PC1*1.1,
           y = p_load$PC2*1.1, 
           label = p_load$rownames,
           colour =  "#494949") +
  labs(caption = str_wrap("Fig. 1: Domesticated African and Asian rice have similar phenotype. 
                    The PC Analysis of the panicle trait dataset reveals that domesticated rice 
                    accession (O. glaberrima and sativa, in green) cluster together separated
                    from wild accession (O. barthii and rufipogon, in blue).
                    The first component explains almost half of the variance and splits wild and
                    cultivated species.
                    This component also splits traits that are intuitively related to high
                    spikelet number, such as branch number and branch length, from traits that 
                    intuitively should correlate inversely with yield, such as
                    internode lenght. The highest loading belongs to secondary branch number and
                    spikelet number, which almost overlap.",
                          width = 60)) +
  xlab(paste0("PC1: ",
              round(var_expl["PC1"], 2),
              "% var.")) +
  ylab(paste0("PC2: ",
              round(var_expl["PC2"], 2),
              "% var.")) +
  theme_bw() +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1),
        legend.text = element_text(face = "italic")) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = p_ch)

pdf(file = "../fig/fig-01-pca-on-pheno.pdf", height = 8)
p_1
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
# 
# summary(lm(pheno_cali$spn ~ pc$x[, 1]))
# ggplot(pheno_cali, aes(x = sbn, y = spn)) +
#   geom_hex(bins = 30) +
#   facet_wrap(facets =  "Species") +
#   theme_bw() +
#   scale_fill_gradient(high = "#001a33", low = "#cce6ff") 
# ggplot(pheno_cali, aes(x = pbn, y = spn)) +
#   geom_hex(bins = 20) +
#   facet_wrap(facets =  "Species")
# ggplot(pheno_cali, aes(x = PanL, y = spn)) +
#   geom_hex(bins = 20) +
#   facet_wrap(facets =  "Species")
# ggplot(pheno_cali, aes(x = PanL, y = pbn)) +
#   geom_hex(bins = 20) +
#   facet_wrap(facets =  "Species")
# ggplot(pheno_cali, aes(x = pbn, y = sbn)) +
#   geom_hex(bins = 20) +
#   facet_wrap(facets =  "Species") +
#   theme_bw() +
#   scale_fill_gradient(high = "#001a33", low = "#cce6ff") 
# 
# pheno_cali %>%
#   split(.$Species) %>%
#   map(~cor(.$pbn, .$sbn))
# 
# 
# pheno_cali %>%
#   split(.$Species) %>%
#   map(~c(range(.$pbn), range(.$sbn)))
# 
# cor(pheno_cali$spn, pheno_cali$sbn)
# 
# plot(pheno_cali$spn ~ pheno_cali$pbn)
# plot(pheno_cali$spn ~ pheno_cali$sbn)
# 
# pheno_mnp <- pheno_mnp %>%
#   filter(complete.cases(.))
# 
# pc_mnp <- pheno_mnp %>%
#   select_at(vars(`Rachis_length (RL)`:`secondary branch internode average length (SbintL in cm)`)) %>%
#   as.data.frame() %>%
#   # filter(complete.cases(.)) %>%
#   prcomp(.,
#          scale. = T,
#          center = T)
# 
# ggplot(as.data.frame(pc_mnp$x)  %>%
#          cbind(Origin = pheno_mnp$Origin),
#        aes(x = PC1, y = PC2, colour = Origin)) +
#   geom_point() +
#   theme_bw()
# 
# pc_mnp$sdev
# 
# plot(pheno_mnp$spn ~ pc_mnp$x[,1])
# plot(pheno_mnp$spn ~ pheno_mnp$sbn)
# plot(pheno_mnp$spn ~ pheno_mnp$pbn)
# 
# pheno_cali %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
#   ggplot(aes(x = Species, y = naxil)) + 
#   geom_boxplot()
# 
# 
# pheno_cali %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
#   ggplot(aes(x = spn, y = naxil, colour = Species)) + 
#   geom_point(alpha = .5) +
#   # facet_wrap(facets = "Species") +
#   theme_bw()
# 
# pheno_cali %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn +TbN)) %>%
#   split(.$Species) %>%
#   map(~mean(.$naxil))
# 
# 
# pheno_mnp %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn + `Tb_nb (TbN)`)) %>%
#   ggplot(aes(x = spn, y = naxil, colour = Species)) + 
#   geom_point(alpha = .5) +
#   # facet_wrap(facets = "Species") +
#   theme_bw()
# 
# pheno_mnp %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn + `Tb_nb (TbN)`)) %>%
#   split(.$Species) %>%
#   map(~mean(.$naxil))
