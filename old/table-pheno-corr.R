library(tidyverse)
source("helper-functions.R")
load("../data/phenotypes.Rdata")
dds <- readRDS("../data-raw/dds.Rds")


# Use median of phenotypic data --------------------------------------------

acc <- c(B88 = "barthii",
         IR64 = "indica",
         Niponbarre = "japonica",
         Tog5681 = "glaberrima",
         W1654 = "rufipogon")

pheno_mnp <- pheno_mnp %>%
  mutate(species = unname(acc[species]))

pheno_mnp <- pheno_mnp %>% 
  group_by(species) %>%
  summarise_at(vars(pbn, sbn, spn), median) %>% 
  ungroup() 


# And the 20.000 genes that pass Deseq2 treshold --------------------------

dat <- get_expression(rownames(dds),
                      dds = dds)
dat <- dat %>%
  group_by(species, locus_id, stage) %>% 
  summarise(median_exp = median(`normalized expression`)) %>%
  spread(key = stage, value = median_exp) %>%
  mutate(fc = (SM + .1) / (PBM + .1),
         log_fc = log2((SM + .1) / (PBM + .1))) %>%
  ungroup() %>%
  mutate(species = as.character(species)) %>%
  # Merge phenotype
  inner_join(pheno_mnp, by = "species") %>%
  arrange(locus_id, species)

# Get correlations (pearson) ----------------------------------------------
cors <- dat %>% 
  group_by(locus_id) %>%
  summarise(cor_fc_pbn = cor(fc, pbn, method = "pearson"),
            cor_fc_sbn = cor(fc, sbn, method = "pearson"),
            cor_fc_spn = cor(fc, spn, method = "pearson"),
            cor_pbm_pbn = cor(PBM, pbn, method = "pearson"),
            cor_pbm_sbn = cor(PBM, sbn, method = "pearson"),
            cor_pbm_spn = cor(PBM, spn, method = "pearson"),
            cor_logfc_pbn = cor(log_fc, pbn, method = "pearson"),
            cor_logfc_sbn = cor(log_fc, sbn, method = "pearson"),
            cor_logfc_spn = cor(log_fc, spn, method = "pearson"))

# save(cors, file = "data/12-2-correlations.Rdata")
cors %>%
  ungroup() %>%
  # dplyr::rename(IDENTIFIER = locus_id) %>%
  # left_join(mapman, by = "IDENTIFIER") %>%
  arrange(cor_logfc_pbn) %>%
  write.csv2("../tables/table-pheno-corr.csv")

svg("../tables/table-pheno-corr.svg")
dev.off()


# get mapman enrichment ---------------------------------------------------

library(fgsea)
load("../data/mapman.Rdata")


mapman_list <- c(#split(mapman$IDENTIFIER, mapman$level1),
  # split(mapman$IDENTIFIER, mapman$level2),
  split(mapman$IDENTIFIER, mapman$level3))
#split(mapman$IDENTIFIER, mapman$level4)

colnames(cors[, -1]) %>% walk(function(i) {
  svg(paste0("../fig/fig-tmp-", 
             str_replace_all(i, "_", "-")
             , ".svg"),
      width = 18,
      height = 10)
  cor_fc_pbn <- test_gsea(setNames(object = cors[, i, drop = TRUE],
                                   nm = cors$locus_id),
                          mapman_list, plot_top = T)
  dev.off()
})

