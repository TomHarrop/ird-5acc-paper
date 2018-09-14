library(tidyverse)
library(fgsea)
dds <- readRDS("../data-raw/dds.Rds")
source("helper-functions.R")


# Load PCA and TF families ------------------------------------------------

# load("../data/all-sep-deseq.Rdata")
load("../data/rlog-pca.Rdata")

# load("../data/func-anno.Rdata") 
# annos <- annos %>% 
#   dplyr::rename(locus_id = "MSU") %>%
#   select(symbol, locus_id)

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

pcx <- pcro %>%
  as.data.frame(.) %>%
  left_join(tf_fam) %>%
  arrange(desc(PC5)) %>%
  mutate(rank_pc5 = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))


# Test Enrichment ---------------------------------------------------------

gsea_pc5 <- test_gsea(set_names(pcx$PC5, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
  dplyr::rename(Family = "pathway")

pcx_tf <- pcx %>%
  left_join(gsea_pc5) %>%
  filter(Family != "none")


# plot enrichement --------------------------------------------------------

p_enr <- pcx_tf %>%
  arrange(padj) %>%
  mutate(facet = paste0(Family,
                        ", adjusted p-value = ",
                        round(padj, 5))) %>%
  mutate(facet = as_factor(facet)) %>%
  select(rank_pc5, PC5, facet) %>%
  View()
  ggplot(aes(x = rank_pc5,
             y = PC5)) +
  geom_linerange(aes(ymin = 0, ymax = PC5), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  facet_grid(facet ~ .) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))

pdf("../fig/suppl-fig-05-tfs-of-pc5.pdf",
    width = 5, height = 50)
print(p_enr)
dev.off()
