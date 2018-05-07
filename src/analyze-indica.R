library(tidyverse)
library(DESeq2)

dds <- readRDS("../data-raw/dds.Rds")

spc <- "indica"

dds <- dds[, dds@colData$accession == spc &  colnames(dds) != "I5"]
cdata <- colData(dds)[dds@colData$accession == spc, -5]

c_indica <- DESeq2::counts(dds, normalized = FALSE)

dds_indica <- DESeqDataSetFromMatrix(countData = c_indica,
                              colData = cdata,
                              design = ~ stage)

dds_indica <- DESeq(dds_indica)

tt_indica <- as.data.frame(results(dds_indica)) %>%
  mutate(locus_id = rownames(.)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

save(tt_indica, dds_indica, file = "../data/indica.Rdata")
                    
# 
# load("../data/func-anno.Rdata")
# 
# annos <- annos %>%
#   dplyr::rename(locus_id = "MSU")
# 
# tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
#   dplyr::rename(locus_id = "Protein.ID") %>%
#   filter(locus_id %in% rownames(dds)) %>%
#   filter(!duplicated(locus_id))
# 
# tst <- tt_indica[tt_indica$padj < .05, ] %>% 
#   left_join(annos) %>%
#   left_join(tf_fam)
