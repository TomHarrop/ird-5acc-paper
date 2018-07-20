library(tidyverse)
library(readxl)
library(pheatmap)

# Function for Geometric Mean ---------------------------------------------

gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Function that reads fluidigm --------------------------------------------

read_fluidigm <- function(path) {
  out <- read_excel(path) %>%
    select(sampleName, GeneName,
           LOC_Name, Value, Type__1) %>%
    rename(sample_name = sampleName,
           target_name = GeneName,
           locus_id = LOC_Name,
           ct_value = Value,
           type = Type__1) %>%
    mutate(ct_value = as.numeric(ct_value)) %>%
    # remove standard curves
    filter(sample_name != "H20",
           !grepl("Mix", sample_name))
  out 
}


# load fluidigms -------------------------------------------------------

fluidigms <- c(alog_exp = "../data-raw/ALOG_CHIP3.xlsx",
               ref_exp = "../data-raw/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)

# Analyze normalizers -----------------------------------------------------

norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))


# merge normalizers and estimate expression -------------------------------

alog_exp <- fluidigms$alog_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))

# Scale -------------------------------------------------------------------

# Store descriptive names for stages

stage_names <- c(N1 = "Rachis_Meristem",
                 N2 = "Primary_Branch_Meristem",
                 N3 = "Spikelet_meristem",
                 N4 = "Young_Flower_Development")

# Store descriptive names for species

species_names <- c(B = "Barthii",
                   G = "Glaberrima",
                   I = "Indica",
                   J = "Japonica",
                   R = "Rufipogon")

# Scale gene expression and scale gene expression separating 

alog_exp <- alog_exp %>%
  # make a variable that encodes the species
  mutate(species = str_sub(sample_name, 1, 1)) %>%
  # make a variable that encodes for the stage
  mutate(stage = str_sub(sample_name, 3, 4)) %>%
  # make a variable that encode replicate
  mutate(replicate = str_sub(sample_name, 6, 6)) %>%
  # scale expression by gene
  group_by(target_name) %>%
  mutate(expr_scale_gene = scale(expression)) %>%
  ungroup() %>%
  # scale expression by geme end species
  group_by(target_name, species) %>%
  mutate(expr_scale_gene_spec = scale(expression)) %>%
  ungroup() %>%
  mutate(expr_scale_gene_spec = ifelse(is.na(expr_scale_gene_spec), 0, expr_scale_gene_spec)) %>%
  # last, save descriptive names for stages and species
  mutate(stage = stage_names[stage]) %>%
  mutate(stage = factor(stage, 
                        levels = c("Rachis_Meristem",
                                   "Primary_Branch_Meristem",
                                   "Spikelet_meristem",
                                   "Young_Flower_Development"))) %>%
  mutate(species = species_names[species])


# Plot heatmap ------------------------------------------------------------

alog_heat <- alog_exp %>%
  # use median of spec-scaled
  group_by(target_name, species, stage) %>%
  summarise(median_expr = median(expr_scale_gene_spec)) %>%
  ungroup() %>% 
  # spread for heatmap
  mutate(species_stage = paste(species, stage, sep = "_")) %>%
  select(target_name, species_stage, median_expr) %>% 
  spread(key = species_stage, value = median_expr) %>%
  # reorder columns
  select(target_name,
         Japonica_Rachis_Meristem,
         Japonica_Primary_Branch_Meristem,
         Japonica_Spikelet_meristem,
         Japonica_Young_Flower_Development,
         Barthii_Rachis_Meristem,
         Barthii_Primary_Branch_Meristem,
         Barthii_Spikelet_meristem,
         Barthii_Young_Flower_Development,
         Glaberrima_Rachis_Meristem,
         Glaberrima_Primary_Branch_Meristem,
         Glaberrima_Spikelet_meristem,
         Glaberrima_Young_Flower_Development,
         Rufipogon_Rachis_Meristem,
         Rufipogon_Primary_Branch_Meristem,
         Rufipogon_Spikelet_meristem,
         Rufipogon_Young_Flower_Development,
         Indica_Rachis_Meristem,
         Indica_Primary_Branch_Meristem,
         Indica_Spikelet_meristem,
         Indica_Young_Flower_Development) %>%
  as.data.frame(.)

pdf("../fig/alog_fluidigm.pdf")
pheatmap(mat = alog_heat %>% select(-target_name),
         # scale = "row",
         # color = colorRampPalette(c("navy", "white", "goldenrod"))(50),
         # color = colorRampPalette(c( "white", "blue4"))(50),
         color = viridis_pal()(50),
         cluster_cols = F,
         gaps_col = 1:5*4,
         cellwidth = 9,
         cellheight = 9,
         labels_row = alog_heat$target_name)
dev.off()
