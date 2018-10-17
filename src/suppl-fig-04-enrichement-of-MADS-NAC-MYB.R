library(tidyverse)
library(magrittr)
library(fgsea)
library(pheatmap)
# library(gridExtra)
library(viridis)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")
filter_abs_pc <- .00185

set.seed(1)


# Load PCA and TF families ------------------------------------------------

# load("../data/all-sep-deseq.Rdata")
load("../data/rlog-pca.Rdata")

load("../data/func-anno.Rdata") 
annos <- annos %>% 
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

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


# Load subfamily definition -----------------------------------------------

# mads
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-8-242

# mads_arora <- read_excel("../data-raw/mads_subfam-copyandpasted_arora_2007.xlsx") %>%
#   select(locus_id, subfamily)

# How many genes pass cutoff ----------------------------------------------



# pcx %>% 
#   # filter(PC5 > filter_abs_pc | PC5 < -filter_abs_pc) %>%
#   filter(PC5 > filter_abs_pc) %>%
#   nrow()

# Test Enrichment ---------------------------------------------------------

gsea_pc5 <- test_gsea(set_names(pcx$PC5, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
  dplyr::rename(Family = "pathway")

pcx_tf <- pcx %>%
  left_join(gsea_pc5) %>%
  filter(Family != "none")



# Save families -----------------------------------------------------------

if(FALSE) {
write_family <- function(family_id) {
  pcx_tf %>%
    filter(Family == family_id) %>%
    select(locus_id, Family) %>%
    mutate(symbol = oryzr::LocToGeneName(locus_id) %$%
             symbols) %>%
    write_csv(paste0("../fix_gene_names/",
                     family_id,
                     ".csv"))
}

write_family("MADS")
write_family("NAC")
write_family("SBP")             

}
  
  pcx_tf %>%
    filter(Family == "MADS") %>%
    select(locus_id, Family) %>%
    mutate(symbol = oryzr::LocToGeneName(locus_id) %$%
             symbols) %>%
    write_csv("../fix_gene_names/MADS.csv")

# Plot Enrichement --------------------------------------------------------

families <- c(mads = "MADS",
              nac = "NAC",
              # myb = "MYB",
              sbp = "SBP")

p_enr <- ggplot(pcx_tf %>%
                  filter(Family %in% families) %>%
                  arrange(padj) %>%
                  mutate(relevant = case_when(PC5 > filter_abs_pc ~ "Yes",
                                              PC5 < -filter_abs_pc ~ "Yes",
                                              TRUE ~ "No"),
                         Family = as_factor(Family)) %>%
                  mutate(facet = paste0(Family,
                                        ", adjusted p-value = ",
                                        round(padj, 4))) %>%
                  mutate(facet = as_factor(facet)),
                aes(x = rank_pc5,
                    y = PC5,
                    colour = relevant)) + 
  geom_linerange(aes(ymin = 0, ymax = PC5), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(limits = c(1, max(pcx_tf$rank_pc5)),
                     breaks = c(1, 5000, 10000, 15000, 20000, max(pcx_tf$rank_pc5)),
                     labels = c("1\n [BM]", "5000", "10000", "15000", "20000",
                                paste0(max(pcx_tf$rank_pc5), "\n [SM]"))) +
  facet_grid(. ~ facet) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6)) +
  labs(x = "Ranks of genes on PC5",
       y = "PC5 Value")

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
    # left_join(subfams) %>% #### ADD SUBFAMILY!
    # left_join(annos) %>% #### ADD GENE NAME
    left_join(pcx_tf %>% select(locus_id, PC5, rank_pc5)) %>%
    as.data.frame() %>%
    # mutate(locus_id = paste(rank_pc5, locus_id, symbol)) %>%
    mutate(locus_id = paste(rank_pc5, locus_id)) %>%
    # select(-symbol) %>%
    arrange(desc(PC5)) %>%
    distinct() %>%
    column_to_rownames("locus_id")
  
  rows_cut <- to_heat %>% 
    mutate(counter = case_when(PC5 > 0 ~ 1,
                               TRUE ~ 0)) %$%
    sum(counter)
  print(rows_cut)
  
  p <- pheatmap(to_heat %>% select(-PC5, -rank_pc5),
                color = viridis_pal()(50),
                show_rownames = T,
                cutree_cols = 2,
                cluster_cols = F,
                cluster_rows = F,
                gaps_col = 5,
                gaps_row = rows_cut,
                cellwidth = 9,
                cellheight = 5,
                main = family,
                fontsize_row = 5)
                # annotation_row = to_heat[, "subfamily", drop = F]) #### ADD SUBFAMILY!
  
  return(p)
}

# Plot families -----------------------------------------------------------

# heat_ap2 <- plot_heatmap("AP2-EREBP")
heat_mads <- plot_heatmap("MADS")
# plot_heatmap("WRKY")
# heat_hb <- plot_heatmap("HB")
heat_nac <- plot_heatmap("NAC")
heat_sbp <- plot_heatmap("SBP", cutree_rows = 1)
# heat_tcp <- plot_heatmap("TCP", cutree_rows = 1)
# heat_bhlh <- plot_heatmap("bHLH")
# heat_zfhd <- plot_heatmap("zf-HD", cutree_rows = 1)

# tst <- families %>% map(plot_heatmap)

# Save plots --------------------------------------------------------------

pdf("../fig/suppl-fig-NAC-MADS-SPL-heatmap.pdf",
    height = 6.2,
    width = 11,
    paper = "a4r")
heats <- cowplot::plot_grid(heat_mads[[4]],
                            heat_nac[[4]],
                            heat_sbp[[4]],
                            labels = c("B", "C", "D"),
                            ncol = 3)
cowplot::plot_grid(p_enr, heats,
                   nrow = 2,
                   rel_heights = c(2,3),
                   labels = c("A", "")) %>%
  cowplot::add_sub(., str_wrap("While MADS and SBP genes are expressed 
                               preferentially in the SM, NAC are more
                               expressed in the BM.
                               A: Enrichment is estimated with the
                               the GSEA method, which returns a permutation based
                               adjusted pvalue of 0.0037 for MADS box genes,
                               of 0.0.0055 for NAC genes and of 0.0245 for SBP genes.
                               B-D: For the heatmap, we used the 10% of genes that
                               have the highest absolute loading on PC5 (shown
                               in red in the enrichment plot).
                               Heatmaps display normalized RNAseq counts
                               which have been z-score scaled independently for each species",
                               width = 80)) %>%
  cowplot::ggdraw() %>% print()

dev.off()




