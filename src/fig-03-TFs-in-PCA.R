library(tidyverse)
library(fgsea)
library(pheatmap)
library(gridExtra)
library(viridis)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")


# Load PCA and TF families ------------------------------------------------

# load("../data/all-sep-deseq.Rdata")
load("../data/rlog-pca.Rdata")

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))


# Prepare data - merge PC and TF ------------------------------------------

# pcx <- pc_spc$x %>%
pcx <- pcro %>%
  as.data.frame(.) %>%
  # rownames_to_column() %>%
  # dplyr::rename(locus_id = "rowname") %>%
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


# Plot Enrichement --------------------------------------------------------

families <- c(ap2 = "AP2-EREBP", 
              mads = "MADS",
              # nac = "NAC",
              hb = "HB")#,
              # sbp = "SBP",
              # bhlh = "bHLH")

p_enr <- ggplot(pcx_tf %>%
                filter(Family %in% families) %>%
                  mutate(facet = paste0(Family,
                                        ", adjusted p-value = ",
                                        round(padj, 3))) %>%
                  mutate(facet = as_factor(facet)),
                aes(x = rank_pc5,
                    y = PC5)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC5), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  # geom_rug(aes(x = rank_pc1, y = NULL), alpha = .5) +
  # facet_grid(Family ~ .) +
  facet_grid(. ~ facet) +
  theme_bw()
# p_enr

# Define functions for heatmap --------------------------------------------

plot_heatmap <- function(family = "AP2-EREBP",
                         norm = by_locus_species,
                         cutree_rows = 2,
                         filter_abs_pc = .002) {
  norm <- enquo(norm)
  to_heat <- pcx_tf %>%
    filter(Family == family) %>%
    mutate(abs_pc5 = abs(PC5)) %>%
    # arrange(desc(abs_pc1)) %>%
    filter(abs_pc5 > filter_abs_pc) %>%
    .$locus_id %>%
    get_expression(dds) %>%
    group_by(locus_id, species, stage) %>%
    summarise(to_plot = median(!!norm)) %>%
    ungroup() %>%
    mutate(stage = as.character(stage)) %>%
    mutate(stage = case_when(stage == "PBM" ~ "BM",
                             TRUE ~ stage)) %>%
    mutate(stage_species = paste(stage, species, sep = "_")) %>%
    select(locus_id, stage_species, to_plot) %>%
    spread(key = stage_species, value = to_plot) %>%
    as.data.frame() %>%
    column_to_rownames("locus_id") 
  
  
  p <- pheatmap(to_heat,
           color = viridis_pal()(50),
           show_rownames = T,
           # fontsize = 5,
           cutree_cols = 2,
           cluster_cols = F,
           # cluster_rows = F,
           gaps_col = 5,
           cutree_rows = cutree_rows,
           cellwidth = 9,
           cellheight = 5,
           main = family,
           fontsize_row = 5)
  
  return(p)
}

# Plot families -----------------------------------------------------------

heat_ap2 <- plot_heatmap("AP2-EREBP")
heat_mads <- plot_heatmap("MADS")
# plot_heatmap("WRKY")
heat_hb <- plot_heatmap("HB")
heat_nac <- plot_heatmap("NAC")
heat_sbp <- plot_heatmap("SBP", cutree_rows = 1)
heat_tcp <- plot_heatmap("TCP", cutree_rows = 1)
heat_bhlh <- plot_heatmap("bHLH")
heat_zfhd <- plot_heatmap("zf-HD", cutree_rows = 1)

# tst <- families %>% map(plot_heatmap)

# Save plots --------------------------------------------------------------

pdf("../fig/fig-05-TFs-in-PCA-locusid.pdf",
    height = 6,
    width = 10)
grid.arrange(grobs = list(p_enr,
                          heat_ap2[[4]],
                          heat_mads[[4]],
                          heat_hb[[4]]),#,
                          # heat_nac[[4]],
                          # heat_sbp[[4]],
                          # heat_bhlh[[4]]),
             # ncol = 4)
             layout_matrix = rbind(c(1,1,1),
                                   c(1,1,1),
                                   2:4,
                                   2:4,
                                   2:4,
                                   2:4,
                                   # 2:4,
                                   2:4))
             # layout_matrix = cbind(c(1,1,1,1),
             #                       c(2,2,5,5),
             #                       c(3,3,6,6),
             #                       c(4,4,7,7)))

dev.off()


# Save Locus_ids ----------------------------------------------------------
# fam_to_csv <- function(fam = "AP2-EREBP") {
#   pcx_tf %>%
#     filter(Family == fam) %>%
#     mutate(abs_pc1 = abs(PC1)) %>%
#     # arrange(desc(abs_pc1)) %>%
#     filter(abs_pc1 > 2) %>%
#     select(locus_id) %>%
#     write_excel_csv(path = paste0("../tables/table-tfs-fig04-",
#                                   fam,
#                                   ".csv"))
# }
# 
# walk(.x = c("AP2-EREBP", "HB", "MADS"),
#      .f = fam_to_csv)


# suppl figure 1 ----------------------------------------------------------

p_enr <- ggplot(pcx_tf %>%
                  arrange(padj) %>%
                  mutate(facet = paste0(Family,
                                        ", adjusted p-value = ",
                                        round(padj, 3))) %>%
                  mutate(facet = as_factor(facet)),
                aes(x = rank_pc5,
                    y = PC5)) +
  geom_linerange(aes(ymin = 0, ymax = PC5), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  facet_wrap(facets = "facet", ncol = 3) +
  theme_bw() 

pdf("../fig/suppl-fig-01-tfs-of-pc5.pdf",
    width = 9, height = 30)
print(p_enr)
dev.off()
