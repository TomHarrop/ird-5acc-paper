library(tidyverse)
library(fgsea)
library(pheatmap)
library(gridExtra)
library(viridis)
source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")
filter_abs_pc <- .002

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
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC5), lwd = 1) + 
  geom_hline(yintercept = 0,
             lwd = .05,
             colour = "grey") +
  scale_color_manual(values = c("black", "red")) +
  # annotate("rect",
  #          fill = "red",
  #          xmin = -Inf, xmax = Inf,
  #          ymin = filter_abs_pc, ymax = Inf,
  #          alpha = .05) +
  # annotate("rect",
  #          fill = "red",
  #          xmin = -Inf, xmax = Inf,
  #          ymin = -Inf, ymax = -filter_abs_pc,
  #          alpha = .05) +
  # geom_rug(aes(x = rank_pc1, y = NULL), alpha = .5) +
  # facet_grid(Family ~ .) +
  facet_grid(. ~ facet) +
  theme_bw()
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

hb_shain <- read_csv("../data-raw/hb_genes_jain2008.csv") %>%
  dplyr::rename(subfamily = "class",
                locus_id = "msuId") %>%
  select(locus_id, subfamily) 

mads_arora <- read_excel("../data-raw/mads_subfam-copyandpasted_arora_2007.xlsx") %>%
  select(locus_id, subfamily)

subfams <- bind_rows(ap2_sharoni, hb_shain, mads_arora)


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
    left_join(subfams) %>% #### ADD SUBFAMILY!
    left_join(annos) %>% #### ADD GENE NAME
    as.data.frame() %>%
    mutate(locus_id = paste(locus_id, symbol)) %>%
    select(-symbol) %>%
    distinct() %>%
    # print(.)
    column_to_rownames("locus_id")
  # 
  # print(data.frame(tst = rownames(to_heat))) %>%
  #   write_csv(path = "~/Desktop/mads.csv")
  
  p <- pheatmap(to_heat %>% select(-subfamily),
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
           fontsize_row = 5,
           annotation_row = to_heat[, "subfamily", drop = F]) #### ADD SUBFAMILY!
  
  return(p)
}

# Plot families -----------------------------------------------------------

heat_ap2 <- plot_heatmap("AP2-EREBP")
heat_mads <- plot_heatmap("MADS")
# plot_heatmap("WRKY")
heat_hb <- plot_heatmap("HB")
# heat_nac <- plot_heatmap("NAC")
# heat_sbp <- plot_heatmap("SBP", cutree_rows = 1)
# heat_tcp <- plot_heatmap("TCP", cutree_rows = 1)
# heat_bhlh <- plot_heatmap("bHLH")
# heat_zfhd <- plot_heatmap("zf-HD", cutree_rows = 1)

# tst <- families %>% map(plot_heatmap)

# Save plots --------------------------------------------------------------

pdf("../fig/fig-05-TFs-in-PCA-locusid_subfams.pdf",
    # height = 6,
    # width = 10)
    height = 6,
    width = 14)
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
