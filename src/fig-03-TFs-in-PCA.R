library(tidyverse)
library(magrittr)
library(fgsea)
library(pheatmap)
library(gridExtra)
library(viridis)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")
filter_abs_pc <- .00185
# filter_abs_pc <- .002

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


# Plot Enrichement --------------------------------------------------------

families <- c(ap2 = "AP2-EREBP", 
              # mads = "MADS",
              # nac = "NAC",
              hb = "HB")#,
              # sbp = "SBP",
              # bhlh = "bHLH")

p_enr <- ggplot(pcx_tf %>%
                  filter(Family %in% families) %>%
                  mutate(relevant = case_when(PC5 > filter_abs_pc ~ "Yes",
                                              PC5 < -filter_abs_pc ~ "Yes",
                                              TRUE ~ "No")) %>%
                  mutate(facet = paste0(Family,
                                        ", adjusted p-value = ",
                                        round(padj, 3))) %>%
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
  labs(x = "Ranks of genes on PC5",
       y = "PC5 Value")

# p_enr


# Load subfamilies --------------------------------------------------------

ap2_sharoni_path <- "../data-raw/ap2s-sharoni/Supplementary_Table_S1.xls"

ap2_sharoni <- bind_cols(read_excel(ap2_sharoni_path,
                                    range = "A4:A166",
                                    col_names = "locus_id"),
                         read_excel(ap2_sharoni_path,
                                    range = "N4:O166",
                                    col_names = c("subgroub", "subfamily"))) %>%
  mutate(locus_id = paste0("LOC_", locus_id),
         subfam = paste(subfamily, subgroub, sep = " - ")) %>%
  select(locus_id, subfamily)
  # left_join(annos %>% select(locus_id, symbol))


ap2_symbols <- ap2_sharoni$locus_id %>%
  oryzr::LocToGeneName() %>%
  filter(!duplicated(MsuID)) %>%
  rename(locus_id = "MsuID",
         symbol = "symbols") %>%
  select(locus_id, symbol)

ap2_sharoni %<>%
  left_join(ap2_symbols)

hb_shain <- read_csv("../data-raw/hb_genes_jain2008.csv") %>%
  dplyr::rename(subfamily = "class",
                locus_id = "msuId") %>%
  select(locus_id, subfamily, symbol) 

# mads_arora <- read_excel("../data-raw/mads_subfam-copyandpasted_arora_2007.xlsx") %>%
#   select(locus_id, subfamily)

subfams <- bind_rows(ap2_sharoni,
                     # mads_arora,
                     hb_shain)

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
    left_join(subfams) %>% #### ADD SUBFAMILY!
    # left_join(annos) %>% #### ADD GENE NAME
    left_join(pcx_tf %>% select(locus_id, PC5, rank_pc5)) %>%
    as.data.frame() %>%
    mutate(locus_id = paste(rank_pc5, locus_id, symbol)) %>%
    select(-symbol) %>%
    arrange(desc(PC5)) %>%
    distinct() %>%
    column_to_rownames("locus_id")

  rows_cut <- to_heat %>% 
    mutate(counter = case_when(PC5 > 0 ~ 1,
                               TRUE ~ 0)) %$%
    sum(counter)
  print(rows_cut)
  
  p <- pheatmap(to_heat %>% select(-subfamily, -PC5, -rank_pc5),
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
           fontsize_row = 5,
           annotation_row = to_heat[, "subfamily", drop = F]) #### ADD SUBFAMILY!
  
  return(p)
}

# Plot families -----------------------------------------------------------

heat_ap2 <- plot_heatmap("AP2-EREBP")
# heat_mads <- plot_heatmap("MADS")
# plot_heatmap("WRKY")
heat_hb <- plot_heatmap("HB")
# heat_nac <- plot_heatmap("NAC")
# heat_sbp <- plot_heatmap("SBP", cutree_rows = 1)
# heat_tcp <- plot_heatmap("TCP", cutree_rows = 1)
# heat_bhlh <- plot_heatmap("bHLH")
# heat_zfhd <- plot_heatmap("zf-HD", cutree_rows = 1)

# tst <- families %>% map(plot_heatmap)

# Save plots --------------------------------------------------------------

pdf("../fig/fig-03-TFs-in-PCA-locusid_subfams.pdf",
    height = 8,
    width = 12)
heats <- cowplot::plot_grid(heat_ap2[[4]],
                            heat_hb[[4]],
                            labels = c("B", "C"))
cowplot::plot_grid(p_enr, heats,
                   nrow = 2,
                   rel_heights = c(2,5),
                   labels = c("A", "")) %>%
  cowplot::add_sub(., str_wrap("AP2-EREBP and HB transcription factor change expression 
                               between BM and SM. 
                               A. AP2 and HB genes are distributed toward the extremes ofPC5 
                               (with adjusted p-values of 0.004 and 0.021 from a GSEA test).
                               We set an experimental cutoff that retains approximately 10% of
                               genes from any of the two extremes (in red).
                               B. Most of AP2-EREBP genes that pass the cutoff are more expressed
                               in the BM. Although only four AP2-EREBP genes are more expressed in
                               the SM, three of them belong to the AP2 subfamily. Instead, genes
                               at the opposite extreme of PC5 belong preferentially to RAV, DREB
                               and ERF subfamilies.
                               C. HB genes that pass the cutoff are prefertially more expressed
                               in the SM. Those genes belong preferentially to the HD-ZIP IV 
                               subfamily",
                   width = 80)) %>%
  cowplot::ggdraw()

dev.off()



