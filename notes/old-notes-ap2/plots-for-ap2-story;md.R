library(tidyverse)
library(DESeq2)
source("../src/helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

"LOC_Os03g05510" %>% 
  get_expression(dds) %>%
  plot_norm_expr()

"LOC_Os09g28440" %>% 
  get_expression(dds) %>%
  plot_norm_expr()

"LOC_Os04g38720"  %>% 
  get_expression(dds) %>%
  plot_norm_expr()

"LOC_Os04g49230" %>% 
  get_expression(dds) %>%
  plot_norm_expr()
  
"LOC_Os02g45850" %>%
  get_expression(dds) %>%
  plot_norm_expr()

