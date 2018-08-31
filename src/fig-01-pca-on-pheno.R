library(tidyverse)
library(ggfortify)
library(grid)
library(gridExtra)
library(ggpubr)  
library(viridis)

load("../data/phenotypes.Rdata")
color_palette <- viridis(20)[c(6, 14, 8, 16)]
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

p_load <- p$layers[[3]]$data

p_1 <- p$data %>%
  ggplot(aes(x = PC1,
             y = PC2,
             fill = Species,
             colour = Species,
             pch = Species)) +
  geom_point(size = 1.2,
             alpha = .7) +
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
print(p_1)
dev.off()

# PCA with labels ----------------------------------------------------------

pc <- pheno_cali %>%
  select_at(vars(Name, RL:spn)) %>%
  # select(-PanL) %>%
  as.data.frame() %>%
  group_by(Name) %>%
  mutate(cnt = 1:n()) %>%
  ungroup() %>%
  mutate(cnt = paste(Name, cnt, sep = "-")) %>%
  as.data.frame() %>%
  column_to_rownames("cnt") %>%
  select(-Name) %>%
  prcomp(.,
         scale. = T,
         center = T)

p <- autoplot(pc, data = pheno_cali, 
              loadings = TRUE,
              loadings.label = TRUE) 

pos <- position_jitter(width = .4,
                       height = 0, 
                       seed = 1)

p$data %>%
  mutate(sequenced = case_when(Name %in% pheno_mnp$species ~ "in rnaseq",
                               TRUE ~ "no")) %>%
  ggplot(aes(x = Origin, 
             y = PC1)) +
  geom_boxplot(outlier.alpha = 0, colour = "black") +
  geom_jitter(data = . %>%
                filter(sequenced == "no"),
              height = 0,
              alpha = .5,
              colour = "darkgrey") +
  geom_point(data = . %>%
               filter(sequenced == "in rnaseq"),
             # alpha = .5,
             colour = "red",
             position = pos) + 
  # scale_color_manual(values = c("red", "grey")) +
  ggrepel::geom_label_repel(data = . %>%
                              filter(sequenced == "in rnaseq"),
                            aes(label = Name),
                            position = pos) +
  theme_bw() 
  
pheno_mnp$species
  

 p$data %>%
  filter(sequenced == "no")