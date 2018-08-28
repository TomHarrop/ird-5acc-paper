library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)

source("helper_functions_fluidigm.R")

# load fluidigms -------------------------------------------------------

fluidigms <- c(alog_exp = "../data-raw/fluidigm_ok/ALOGok_CHIP3.xlsx",
               ap2_exp = "../data-raw/fluidigm_ok/AP2ok_CHIP3.xlsx",
               horm_exp = "../data-raw/fluidigm_ok/HormoneOK_CHIP3.xlsx",
               cluster_exp = "../data-raw/fluidigm_ok/VariousClusterOK_CHIP3.xlsx",
               various_exp = "../data-raw/fluidigm_ok/VariousKnownOK_CHIP3.xlsx",
               ref_exp = "../data-raw/fluidigm/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

# Merge results -----------------------------------------------------------

chip3_exp <- fluidigms; chip3_exp$ref_exp <- NULL
chip3_exp <- chip3_exp %>% purrr::reduce(bind_rows)

# merge normalizers and estimate expression -------------------------------

chip3_exp <- chip3_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))


# Check distribution ------------------------------------------------------

# ggplot(chip3_exp, aes(x = sample_name,
#                       y = expr_scale_gene,
#                       fill = species)) +
#   geom_boxplot() +
#   ggplot2::scale_y_log10()

# Plot both RNAseq and fluidigm -------------------------------------------

library(DESeq2)
library(gridExtra)
library(grid)

source("../src/helper-functions.R")

dds <- readRDS("../data-raw/dds.Rds")

chip3_rnaseq <- unique(chip3_exp$locus_id) %>%
  get_expression(dds = dds) %>%
  mutate(species = factor(species,
                          levels = c("japonica", 
                                     "barthii",
                                     "glaberrima", 
                                     "rufipogon", 
                                     "indica"))) 

chip3_exp <- chip3_exp %>% 
  scale_tidy_fluidigm() %>%
  mutate(species = factor(species,
                          levels = c("Osj", 
                                     "Ob",
                                     "Og", 
                                     "Or", 
                                     "Osi"))) 



# id_fluidigm <- "LOC_Os06g47150"
# rnaseq_dat <- clust_rnaseq

# plot_both("LOC_Os04g36054",
#           fluidigm_dat = clust_exp,
#           rnaseq_dat = clust_rnaseq)

pdf(file = "../fig/chip3-fluidigm-rnaseq.pdf",
    height = 7)
unique(chip3_rnaseq$locus_id) %>% map(~plot_both(.,
                                              fluidigm_dat = chip3_exp,
                                              rnaseq_dat = chip3_rnaseq))
dev.off()


# Load cluster data -------------------------------------------------------

cls <- read_csv(file = "../data-raw/annotated_clusters_scaled_l2fc.csv") %>%
  dplyr::rename(locus_id2 = "MsuID") %>%
  select(locus_id2, cluster)



# Merge cluster and fluidigm ----------------------------------------------

clust_heat <- chip3_exp %>%
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

pdf("../fig/fluidigm_chip3_clusters.pdf",
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


# Interesting ones --------------------------------------------------------

ints <- c(g1l2 = "LOC_Os06g46030",
          rav2 = "LOC_Os01g04800",
          ap37 = "LOC_Os01g58420",
          dreb1g = "LOC_Os02g45450",
          erebp153 = "LOC_Os09g25600",
          erebp86 = "LOC_Os03g19900",
          erf64 = "LOC_Os03g08500",
          erf109 = "LOC_Os09g13940",
          erf142 = "LOC_Os05g32270",
          erf20 = "LOC_Os02g45420",
          erf36 = "LOC_Os10g41130",
          erf47 = "LOC_Os03g09170",
          erf76 = "LOC_Os04g57340",
          erf91 = "LOC_Os02g43790",
          erf130 = "LOC_Os05g41760",
          
          PLATZ = "LOC_Os04g50120",
          APO1 = "LOC_Os06g45460",
          APO2 = "LOC_Os04g51000",
          valid11 = "LOC_Os05g01940",
          valid06 = "LOC_Os05g03760")
