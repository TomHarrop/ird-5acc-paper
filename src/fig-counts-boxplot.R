library(tidyverse)
library(DESeq2)
dds <- readRDS("../data-raw/dds.Rds")
load("../data/rlog-pca.Rdata")


svg("../fig/fig-counts-boxplot.svg",
    height = 8,
    width = 8)
par(mfrow = c(3,1))
boxplot(DESeq2::counts(dds),
        main = "Raw Counts",
        outline = FALSE)
boxplot(DESeq2::counts(dds, normalize = T),
        main = "Normalized Counts",
        outline = FALSE)
boxplot(assay(rld),
        main = "rlog Counts",
        outline = FALSE)
dev.off()
