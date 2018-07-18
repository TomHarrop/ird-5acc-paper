library(tidyverse)
library(ggfortify)
library(gridExtra)

load("../data/phenotypes.Rdata")

# Test and plot three correlations ----------------------------------------

pbn_spn <- pheno_cali %>% 
  ggplot(aes(x = pbn, y = spn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 

sbn_spn <- pheno_cali %>% 
  ggplot(aes(x = sbn, y = spn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 

pbn_sbn <- pheno_cali %>% 
  ggplot(aes(x = pbn, y = sbn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  labs(caption = str_wrap("Correlation between the main traits that
                          determine panicle phenotype. A. Primary Branch 
                          Number and Spikelet Number correlates only slightly
                          mainly in wild species. B. Secondary Branch Number 
                          and Spikelet number highly correlates in cultivated
                          species, less in wild species. C. Primary and
                          Secondary Branches poorly correlates, suggesting that 
                          they are controlled by different genetic mechanisms",
                          width = 70))+
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1))


# Arrange and save plots --------------------------------------------------------------

pdf(file = "../fig/fig-02-descript-pheno-stat.pdf", 
    width = 8,
    height = 9)
grid.arrange(pbn_spn,
             sbn_spn,
             pbn_sbn,
             layout_matrix = matrix(c(1,1 ,2,2,3,3,3), ncol = 1))
dev.off()


# Correlation estimates ---------------------------------------------------

pheno_cali %>%
  group_by(Species) %>%
  summarize(corr_spn_sbn = cor(spn, sbn),
            corr_spn_pbn = cor(spn, pbn),
            corr_pbn_sbn = cor(pbn, sbn))

