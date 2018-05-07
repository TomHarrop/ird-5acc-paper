library(tidyverse)
library(DESeq2)

dds <- readRDS("../data-raw/dds.Rds")


analyze_spc <- function(spc) {
  dds <- dds[, dds@colData$accession == spc]
  cdata <- colData(dds)[dds@colData$accession == spc, -5]
  
  cnts <- DESeq2::counts(dds, normalized = FALSE)
  
  dds <- DESeqDataSetFromMatrix(countData = cnts,
                                       colData = cdata,
                                       design = ~ stage)
  
  dds <- DESeq(dds)
  return(dds)
}

get_signif  <- function(dds)
  {
  res <- as.data.frame(results(dds)) %>%
    mutate(locus_id = rownames(.)) %>%
    filter(!is.na(padj)) %>%
    arrange(padj)
  return(res)
}
  

spc <- unique(as.character(colData(dds)$accession))
spc <- set_names(x = spc, nm = spc)
dds_spc <- spc %>% map(analyze_spc)
spc_res <- dds_spc %>% map(get_signif)

spc_res %>% map(~nrow(filter(., padj < .05)))

save(dds_spc, spc_res, file = "../data/all-sep-deseq.Rdata")

