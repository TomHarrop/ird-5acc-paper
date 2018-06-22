library(tidyverse)
library(fgsea)
source("helper-functions.R")
# dds <- readRDS("../data-raw/dds.Rds")
load("../data/rlog-pca.Rdata")

fams <- read_delim("../data-raw/famInfo.table.txt",
                   delim = "\t") %>%
  select(MSU, Name) %>%
  dplyr::rename(locus_id = "MSU",
                Family = "Name") %>%
  filter(locus_id != "none")


plot_tfs_pc <- function(pc = PC5,
                        pdf_name = "suppl-fig-enrich-fams-funricegene") 
{
  sel_pc <- enquo(pc)
  
  # select pc
  pcx <- pcro %>% 
    left_join(fams) %>%
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
  
  pdf(paste0("../fig/",
             pdf_name,
             deparse(substitute(sel_pc)),
             ".pdf"),
      width = 14, height = 200)
  print(p_enr)
  dev.off()
}

plot_tfs_pc()
