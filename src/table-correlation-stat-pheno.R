load("../data/all-sep-deseq.Rdata")
load("../data/phenotypes.Rdata")
library(pheatmap)
library(tidyverse)

source("helper-functions.R")
# dds <- readRDS("../data-raw/dds.Rds")


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

# dat <- get_expression(rownames(dds),
#                       dds = dds)
# dat <- dat %>%
#   group_by(species, locus_id, stage) %>% 
#   summarise(median_exp = median(`normalized expression`)) %>%
#   spread(key = stage, value = median_exp) %>%
#   mutate(fc = (SM + .1) / (PBM + .1),
#          log_fc = log2((SM + .1) / (PBM + .1))) %>%
#   ungroup() %>%
#   mutate(species = as.character(species)) %>%
#   # Merge phenotype
#   inner_join(pheno_mnp, by = "species") %>%
#   arrange(locus_id, species)

dat <- spc_res_df %>% 
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  gather(stat_barthii:stat_japonica, key = species, value = stat) %>%
  mutate(species = str_split_fixed(string = species, pattern = "_", n = 2)[,2]) %>%
  inner_join(pheno_mnp) %>%
  arrange(locus_id, species)



# Get correlations (pearson) ----------------------------------------------
cors <- dat %>% 
  group_by(locus_id) %>%
  summarise(cor_stat_pbn = cor(stat, pbn, method = "pearson"),
            cor_stat_sbn = cor(stat, sbn, method = "pearson"),
            cor_stat_spn = cor(stat, spn, method = "pearson"))


# get mapman enrichment ---------------------------------------------------

library(fgsea)


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

cors %>% ungroup()


# save(cors, file = "data/12-2-correlations.Rdata")
load("../data/mapman.Rdata")
load("../data/func-anno.Rdata")
dds <- readRDS("../data-raw/dds.Rds")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

mapman <- mapman %>%
  dplyr::rename(locus_id = "IDENTIFIER") %>%
  filter(!duplicated(locus_id))

cors <- cors %>%
  ungroup() %>%
  # dplyr::rename(IDENTIFIER = locus_id) %>%
  # left_join(mapman, by = "IDENTIFIER") %>%
  arrange(cor_stat_sbn) %>%
  left_join(annos) %>%
  left_join(mapman) %>%
  left_join(tf_fam) 

cors %>% write.csv2("../tables/table-correlation-stat-pheno.csv")

pdf('../fig/fig-correlation-stat-sbn-top100.pdf',
    width = 20, height = 100)
p <- cors$locus_id[1:100] %>%
  get_expression(dds) %>%
  left_join(annos) %>%
  left_join(tf_fam) %>%
  left_join(mapman) %>%
  mutate(locus_id = as_factor(locus_id)) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id",
                        "symbol",
                        "Family",
                        "DESCRIPTION"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 50,
                                       multi_line = T))
print(p)
dev.off()
