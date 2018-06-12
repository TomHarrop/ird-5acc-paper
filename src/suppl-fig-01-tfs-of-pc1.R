library(tidyverse)
library(fgsea)
library(pheatmap)
library(gridExtra)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")


# Load PCA and TF families ------------------------------------------------

load("../data/all-sep-deseq.Rdata")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

pcx <- pc_spc$x %>% 
  as.data.frame(.) %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam) %>%
  arrange(desc(PC1)) %>%
  mutate(rank_pc1 = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))

# Test Enrichment ---------------------------------------------------------

gsea_pc1 <- test_gsea(set_names(pcx$PC1, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
  dplyr::rename(Family = "pathway")

pcx_tf <- pcx %>%
  left_join(gsea_pc1) %>%
  filter(Family != "none")


# Plot enrichement --------------------------------------------------------

p_enr <- ggplot(pcx_tf %>%
                  arrange(padj) %>%
                  mutate(facet = paste0(Family,
                                        ", adjusted p-value = ",
                                        round(padj, 3))) %>%
                  mutate(facet = as_factor(facet)),
                aes(x = rank_pc1,
                    y = PC1)) +
  geom_linerange(aes(ymin = 0, ymax = PC1), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  facet_wrap(facets = "facet", ncol = 3) +
  theme_bw() 

pdf("../fig/suppl-fig-01-tfs-of-pc1.pdf",
    width = 9, height = 30)
p_enr
dev.off()
