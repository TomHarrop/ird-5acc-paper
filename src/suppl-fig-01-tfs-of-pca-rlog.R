library(tidyverse)
library(fgsea)
library(pheatmap)
library(gridExtra)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")


# Load PCA and TF families ------------------------------------------------

load("../data/rlog-pca.Rdata")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

plot_tfs_pc <- function(pc = PC5) 
{
  sel_pc <- enquo(pc)
  
  # select pc
  pcx <- pcro %>% 
    left_join(tf_fam) %>%
    arrange(desc(!!sel_pc)) %>%
    mutate(rank_pc = 1:nrow(.)) %>%
    mutate(Family = ifelse(is.na(Family), "none", Family))
  
  # Test Enrichment
  
  gsea_pc <- test_gsea(set_names(pcx %>%
                                   select(!!sel_pc) %>%
                                   simplify() %>%
                                   unname(), 
                                  nm = pcx$locus_id),
                        mapman = split(pcx$locus_id, pcx$Family)) %>%
    dplyr::rename(Family = "pathway")
  
  pcx_tf <- pcx %>%
    left_join(gsea_pc) %>%
    filter(Family != "none")
  
  p_enr <- ggplot(pcx_tf %>%
                    arrange(padj) %>%
                    mutate(facet = paste0(Family,
                                          ", adjusted p-value = ",
                                          round(padj, 3))) %>%
                    mutate(facet = as_factor(facet)),
                  aes(x = rank_pc,
                      y = !!sel_pc)) +
    geom_linerange(aes(ymin = 0, ymax = !!sel_pc), lwd = 1) + 
    geom_hline(yintercept = 0,
               lwd = .05,
               colour = "grey") +
    facet_wrap(facets = "facet", ncol = 3) +
    theme_bw() 

  pdf(paste0("../fig/suppl-fig-01-rlog-tfs-of-",
             deparse(substitute(sel_pc)),
            ".pdf"),
      width = 9, height = 30)
  print(p_enr)
  dev.off()
}

# plot_tfs_pc(PC1)
# plot_tfs_pc(PC2)
# plot_tfs_pc(PC3)
# plot_tfs_pc(PC4)
# plot_tfs_pc()

plot_tfs_pc()

# Plot enrichement --------------------------------------------------------

# p_enr <- ggplot(pcx_tf %>%
#                   arrange(padj) %>%
#                   mutate(facet = paste0(Family,
#                                         ", adjusted p-value = ",
#                                         round(padj, 3))) %>%
#                   mutate(facet = as_factor(facet)),
#                 aes(x = rank_pc5,
#                     y = PC5)) +
#   geom_linerange(aes(ymin = 0, ymax = PC1), lwd = 1) + 
#   geom_hline(yintercept = 0,
#              lwd = .05,
#              colour = "grey") +
#   facet_wrap(facets = "facet", ncol = 3) +
#   theme_bw() 
# 
# pdf("../fig/suppl-fig-01-tfs-of-pc1.pdf",
#     width = 9, height = 30)
# print(p_enr)
# dev.off()
