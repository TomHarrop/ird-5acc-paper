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


# Plot Enrichement --------------------------------------------------------

families <- c("AP2-EREBP", "WRKY", 
              "MADS", "HB")

p_enr <- ggplot(pcx_tf %>%
                filter(Family %in% families),
                aes(x = rank_pc1,
                    y = PC1)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC1), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  # geom_rug(aes(x = rank_pc1, y = NULL), alpha = .5) +
  facet_wrap(facets = "Family", ncol = 1) +
  theme_bw() 


# Define functions for heatmap --------------------------------------------

plot_heatmap <- function(family = "AP2-EREBP",
                         norm = by_locus_species) {
  norm <- enquo(norm)
  to_heat <- pcx_tf %>%
    filter(Family == family) %>%
    mutate(abs_pc1 = abs(PC1)) %>%
    # arrange(desc(abs_pc1)) %>%
    filter(abs_pc1 > 2) %>%
    .$locus_id %>%
    get_expression(dds) %>%
    group_by(locus_id, species, stage) %>%
    summarise(to_plot = median(!!norm)) %>%
    ungroup() %>%
    mutate(stage_species = paste(stage, species, sep = "_")) %>%
    select(locus_id, stage_species, to_plot) %>%
    spread(key = stage_species, value = to_plot) %>%
    as.data.frame() %>%
    column_to_rownames("locus_id") 
  
  
  p <- pheatmap(to_heat,
           color = colorRampPalette(c( "white", "blue4"))(50),
           show_rownames = F,
           cutree_cols = 2,
           cluster_cols = F,
           # cluster_rows = F,
           gaps_col = 5,
           cutree_rows = 2,
           cellwidth = 9,
           cellheight = 5,
           main = family)
  
  return(p)
}

# Plot families -----------------------------------------------------------

heat_ap2 <- plot_heatmap("AP2-EREBP")
heat_mads <- plot_heatmap("MADS")
plot_heatmap("WRKY")
heat_hb <- plot_heatmap("HB")

# Save plots --------------------------------------------------------------

pdf("../fig/fig-04-TFs-in-PCA.pdf",
    height = length(families)*1.2,
    width = 12)
grid.arrange(grobs = list(p_enr,
                          heat_ap2[[4]],
                          heat_hb[[4]],
                          heat_mads[[4]]),
             nrow = 1)

dev.off()


# Save Locus_ids ----------------------------------------------------------
fam_to_csv <- function(fam = "AP2-EREBP") {
  pcx_tf %>%
    filter(Family == fam) %>%
    mutate(abs_pc1 = abs(PC1)) %>%
    # arrange(desc(abs_pc1)) %>%
    filter(abs_pc1 > 2) %>%
    select(locus_id) %>%
    write_excel_csv(path = paste0("../tables/table-tfs-fig04-",
                                  fam,
                                  ".csv"))
}

walk(.x = c("AP2-EREBP", "HB", "MADS"),
     .f = fam_to_csv)
