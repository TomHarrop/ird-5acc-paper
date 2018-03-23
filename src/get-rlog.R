library(DESeq2)

dds <- readRDS("../data/dds.Rds")

rld <- rlog(dds, blind = FALSE)
save(rld, file = "../data/rld.Rdata")