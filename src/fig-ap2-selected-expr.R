library(tidyverse)
library(DESeq2)
source("helper-functions.R")


genes <- c(LOC_Os08g31580 ="ERF48",
           LOC_Os05g32270 = "ERF142",
           LOC_Os07g03250 = "PLT8")

dds <- readRDS("../data/dds.Rds")

expr <- get_expression(locus_ids = names(genes), dds = dds)

expr <- expr %>%
  mutate(locus_id = paste(locus_id,
                          genes[locus_id],
                          sep = "\n"))

svg("../fig/fig-ap2-selected-expr.svg",
        width = 7,
        height = 3.5)
plot_norm_expr(expr)
dev.off()
