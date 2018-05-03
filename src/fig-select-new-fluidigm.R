library(tidyverse)
library(readxl)
library(DESeq2)
source("helper-functions.R")


# selected <- read_excel("../data-raw/candidate_new_fluidigm.xlsx") %>%
#   arrange(locus_id)
selected <- read_csv("../selected_genes/candidate_new_fluidigm.csv") %>%
  arrange(locus_id) %>%
  select(-X1)

dds <- readRDS("../data-raw/dds.Rds")

selected_expr <- get_expression(locus_ids = selected$locus_id,
                           dds = dds) %>%
  left_join(selected)



pdf(file = "../fig/fig-select-new-fluidigm.pdf",
    width = 16, height = 14)
plot_norm_expr(selected_expr) +
  facet_wrap(facets = c("locus_id",
                        "gene_name",
                        "why_interesting"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 42))
dev.off()

