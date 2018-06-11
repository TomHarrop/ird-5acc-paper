library(tidyverse)
source("helper-functions.R")
library(fgsea)


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


# Plot Enrichement --------------------------------------------------------

p_enr <- ggplot(pcx_tf,
                # filter(Family %in% c("AP2-EREBP", "MADS", "WRKY")),
                aes(x = rank_pc1,
                    y = PC1)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC1), lwd = .5) + 
  geom_hline(yintercept = 0, lwd = .05) +
  geom_rug(aes(x = rank_pc1, y = NULL), alpha = .5) +
  facet_wrap(facets = c( "padj", "Family"), ncol = 1) +
  theme_bw()

pdf("../fig/fig-04-TFs-in-PCA.pdf",
    height = 50, width = 6)
print(p_enr)
dev.off()