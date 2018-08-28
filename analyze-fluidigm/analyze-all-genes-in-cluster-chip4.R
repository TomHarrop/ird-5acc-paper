library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")

# load fluidigms -------------------------------------------------------

fluidigms <- c(hb = "../data-raw/fluidigm_ok/HBok_CHIP4.xlsx",
               horm_exp = "../data-raw/fluidigm_ok/HormoneOK_CHIP4.xlsx",
               mads = "../data-raw/fluidigm_ok/MADSok_CHIP4.xlsx",
               spl = "../data-raw/fluidigm_ok/SPok_CHIP4.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP4.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

# Merge results -----------------------------------------------------------

chip4_exp <- fluidigms; chip4_exp$ref_exp <- NULL
chip4_exp <- chip4_exp %>% purrr::reduce(bind_rows)

# merge normalizers and estimate expression -------------------------------

chip4_exp <- chip4_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Plot both RNAseq and fluidigm -------------------------------------------

library(DESeq2)
library(gridExtra)
library(grid)

source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

chip4_rnaseq <- unique(chip4_exp$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "glaberrima", 
                                     "rufipogon", 
                                     "indica"))) 

chip4_exp <- chip4_exp %>% 
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj", 
                                     "Ob",
                                     "Og", 
                                     "Or", 
                                     "Osi"))) 
# missing locus ids

# LOC_Os10g41231, LOC_Os10g41232, LOC_Os10g41233, LOC_Os10g41234,
# LOC_Os10g41235, LOC_Os10g41236, LOC_Os10g41237, LOC_Os10g41238,
# LOC_Os10g41239, LOC_Os10g41241, LOC_Os10g41242, LOC_Os10g41243,
# LOC_Os10g41244, LOC_Os10g41245, LOC_Os10g41246, LOC_Os10g41247,
# LOC_Os10g41248, LOC_Os10g41249, LOC_Os10g41251, LOC_Os10g41252,
# LOC_Os10g41253, LOC_Os10g41254, LOC_Os10g41255, LOC_Os10g41256,
# LOC_Os10g41257, LOC_Os10g41258, LOC_Os10g41259, LOC_Os10g41261,
# LOC_Os10g41262, LOC_Os10g41263, LOC_Os10g41264, LOC_Os10g41265,
# LOC_Os10g41266, LOC_Os10g41267, LOC_Os10g41268, LOC_Os10g41269,
# LOC_Os10g41270, LOC_Os10g41271, LOC_Os10g41272, LOC_Os10g41273,
# LOC_Os10g41274, LOC_Os10g41275, LOC_Os10g41276, LOC_Os10g41277,
# LOC_Os10g41278, LOC_Os10g41279, LOC_Os10g41281, LOC_Os10g41282,
# LOC_Os10g41283, LOC_Os10g41284, LOC_Os10g41285, LOC_Os10g41286, 
# LOC_Os10g41287, LOC_Os10g41288, LOC_Os10g41289, LOC_Os12g01120, 
# LOC_OS02g52340


# Check distribution ------------------------------------------------------

ggplot(chip4_exp, aes(x = sample_name,
                      y = expression,
                      fill = species)) +
  geom_boxplot() +
  ggplot2::scale_y_log10()

# id_fluidigm <- "LOC_Os06g47150"
# rnaseq_dat <- clust_rnaseq


# plot all ----------------------------------------------------------------


pdf(file = "../fig/chip4-fluidigm-rnaseq.pdf",
    height = 7)
unique(chip4_rnaseq$locus_id) %>% map(~plot_both(.,
                                                 fluidigm_dat = chip4_exp,
                                                 rnaseq_dat = chip4_rnaseq))
dev.off()


# Interesting ones --------------------------------------------------------

# Load cluster data -------------------------------------------------------

cls <- read_csv(file = "../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  dplyr::rename(locus_id2 = "MsuID") %>%
  select(locus_id2, cluster)



# Merge cluster and fluidigm ----------------------------------------------

clust_heat <- chip4_exp %>%
  scale_tidy_fluidigm() %>%
  prepare_for_heat_spec() %>%
  mutate(locus_id2 = str_split_fixed(string = locus_id,
                                     pattern = " ",
                                     2)[, 1]) %>%
  inner_join(cls) %>%
  dplyr::arrange(cluster) %>%
  mutate(cluster = as.factor(cluster)) %>%
  column_to_rownames("locus_id2")


# heatmap -----------------------------------------------------------------

get_gaps <- function(iter) {
  out <- numeric()
  for(i in 1:(length(iter)-1)) {
    if(iter[i] != iter[i+1]) {
      out <- c(out, i) 
    }
  }
  return(out)
}

pdf("../fig/fluidigm_chip4_clusters.pdf",
    width = 12)
pheatmap(mat = clust_heat %>%
           select(-locus_id, -cluster),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         cluster_rows = F,
         gaps_col = 1:5*4,
         gaps_row = get_gaps(clust_heat$cluster),
         cellwidth = 9,
         cellheight = 9, 
         annotation_row = clust_heat[, "cluster", drop = FALSE],
         labels_row = clust_heat$locus_id)
dev.off()

