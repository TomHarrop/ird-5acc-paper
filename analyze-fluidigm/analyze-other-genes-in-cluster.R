library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")


# load fluidigms -------------------------------------------------------

fluidigms4 <- c(mads_exp = "../data-raw/fluidigm/MADS_CHIP4.xlsx",
               spl_exp = "../data-raw/fluidigm/SPL_CHIP4.xlsx",
               hb_exp = "../data-raw/fluidigm/HB_CHIP4.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP4.xlsx") %>%
  map(., read_fluidigm)

fluidigms4$gene_exp <- bind_rows(fluidigms4$mads_exp,
                                 fluidigms4$spl_exp,
                                 fluidigms4$hb_exp)

fluidigms3 <- c(ap2_exp = "../data-raw/fluidigm/AP2_CHIP3.xlsx",
                ref_exp = "../data-raw/fluidigm/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms3 <- fluidigms3$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

norms4 <- fluidigms4$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

# merge normalizers and estimate expression -------------------------------

chip3_exp <- fluidigms3$ap2_exp %>%
  left_join(norms3)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))


chip4_exp <- fluidigms4$gene_exp %>%
  left_join(norms4)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))



# MERGE and SELECT-------------------------------------------------------------------

keep <- read_excel("../data-raw/fluidigm/GENEs in cluster.xlsx")

clust_exp <- bind_rows(chip3_exp, chip4_exp) %>%
  filter(locus_id %in% keep$LOC_Name)

# Plot both RNAseq and fluidigm -------------------------------------------

library(DESeq2)
library(gridExtra)
library(grid)

source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

clust_rnaseq <- unique(clust_exp$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "rufipogon", 
                                     "glaberrima", 
                                     "indica"))) 

clust_exp <- clust_exp %>% 
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj", 
                                     "Ob",
                                     "Or", 
                                     "Og", 
                                     "Osi"))) 



# id_fluidigm <- "LOC_Os06g47150"
# rnaseq_dat <- clust_rnaseq

# plot_both("LOC_Os04g36054",
#           fluidigm_dat = clust_exp,
#           rnaseq_dat = clust_rnaseq)

pdf(file = "../fig/cluster-fluidigm-rnaseq-others.pdf",
    height = 7)
unique(clust_exp$locus_id) %>% map(~plot_both(.,
                                              fluidigm_dat = clust_exp,
                                              rnaseq_dat = clust_rnaseq))
dev.off()


# fluidigm heatmap --------------------------------------------------------

clust_exp <- bind_rows(chip3_exp, chip4_exp) %>%
  filter(locus_id %in% keep$LOC_Name) %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat_spec()

keep <- keep %>%
  dplyr::rename(locus_id2 = "LOC_Name")

clust_heat <- clust_exp %>%
  mutate(locus_id2 = str_split_fixed(string = locus_id, pattern = " ", 2)[, 1]) %>%
  left_join(keep %>% select(locus_id2, Cluster)) %>%
  dplyr::arrange(Cluster) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  column_to_rownames("locus_id2")

get_gaps <- function(iter) {
  out <- numeric()
  for(i in 1:(length(iter)-1)) {
    if(iter[i] != iter[i+1]) {
      out <- c(out, i) 
    }
  }
  return(out)
}

pdf("../fig/cluster_fluidigm_others.pdf",
    width = 12)
pheatmap(mat = clust_heat %>%
           select(-locus_id, -Cluster),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         cluster_rows = F,
         gaps_col = 1:5*4,
         gaps_row = get_gaps(clust_heat$Cluster),
         cellwidth = 9,
         cellheight = 9, 
         annotation_row = clust_heat[, "Cluster", drop = FALSE],
         labels_row = clust_heat$locus_id)
dev.off()
