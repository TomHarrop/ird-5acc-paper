library(tidyverse)
library(DESeq2)
library(magrittr)

source("../src/helper-functions.R")
dds <- readRDS(file = "../data-raw/dds.Rds")

cl <- read_csv("../data-raw/annotated_clusters_scaled_l2fc.csv")
  
cl %>%
  mutate(id = paste(RapID, symbols)) %>%
  select(id, cluster, rufipogon, indica, barthii, glaberrima) %>%
  gather(rufipogon:glaberrima, key = species, value = lfc) %>%
  ggplot(aes(x = species, y = lfc)) +
  geom_boxplot() +
  facet_wrap("cluster", nrow = 3) + 
  theme_bw()



cl %>% filter(grepl("AP2", MsuAnnotation)) %>%
  View()
cl %>% filter(grepl("AP2", MsuAnnotation)) %$%
  get_expression(MsuID, dds) %>%
  plot_norm_expr()
