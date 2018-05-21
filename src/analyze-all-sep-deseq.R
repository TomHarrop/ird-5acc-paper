library(tidyverse)
library(DESeq2)

dds <- readRDS("../data-raw/dds.Rds")


# Define functions --------------------------------------------------------

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


# Run deseq on every species ----------------------------------------------


spc <- unique(as.character(colData(dds)$accession))
spc <- set_names(x = spc, nm = spc)
dds_spc <- spc %>% map(analyze_spc)


# Get significant DE for each species --------------------------------------

spc_res <- dds_spc %>% map(get_signif)
spc_res %>% map(~nrow(filter(., padj < .05)))


# PCA on DESeq2 stat ------------------------------------------------------


spc_res_df <- spc_res %>%
  map(~.[, c("locus_id", "stat")]) %>%
  map(as.data.frame) 

change_name <- function(i) {
  colnames(spc_res_df[[i]])[2] <- paste("stat", i, sep = "_")
  return(spc_res_df[[i]])
}

spc_res_df <- names(spc_res_df) %>% map(change_name) 

spc_res_df <- spc_res_df %>% purrr::reduce(full_join, by = "locus_id")

rownames(spc_res_df) <- spc_res_df$locus_id
spc_res_df$locus_id <- NULL
spc_res_df[is.na(spc_res_df)] <- 0

pc_spc <- prcomp(spc_res_df, scale. = T)


# Save all ----------------------------------------------------------------



save(dds_spc, spc_res, pc_spc, file = "../data/all-sep-deseq.Rdata")

