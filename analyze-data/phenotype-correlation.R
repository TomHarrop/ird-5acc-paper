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
# spc_res %>% map(~nrow(filter(., padj < .05)))


# PCA on DESeq2 stat ------------------------------------------------------


spc_res_df <- spc_res %>%
  map(~.[, c("locus_id", "log2FoldChange")]) %>%
  map(as.data.frame) 

change_name <- function(i) {
  colnames(spc_res_df[[i]])[2] <- paste("l2fc", i, sep = "_")
  return(spc_res_df[[i]])
}

spc_res_df <- names(spc_res_df) %>% map(change_name) 

spc_res_df <- spc_res_df %>% purrr::reduce(full_join, by = "locus_id")


# load phenotypes ---------------------------------------------------------

load("../data/phenotypes.Rdata")

# Use median of phenotypic data --------------------------------------------

acc <- c(B88 = "barthii",
         IR64 = "indica",
         Niponbarre = "japonica",
         Tog5681 = "glaberrima",
         W1654 = "rufipogon")

pheno_mnp <- pheno_mnp %>%
  mutate(species = unname(acc[species]))

pheno_mnp <- pheno_mnp %>% 
  group_by(species) %>%
  summarise_at(vars(pbn, sbn, spn), median) %>% 
  ungroup() 


# Estimate correlation ----------------------------------------------------

dat <- spc_res_df %>% 
  # rownames_to_column() %>%
  # dplyr::rename(locus_id = "rowname") %>%
  gather(l2fc_barthii:l2fc_japonica, key = species, value = l2fc) %>%
  mutate(species = str_split_fixed(string = species, pattern = "_", n = 2)[,2]) %>%
  inner_join(pheno_mnp) %>%
  arrange(locus_id, species)

cors <- dat %>% 
  group_by(locus_id) %>%
  summarise(cor_stat_pbn = cor(l2fc, pbn, method = "pearson"),
            cor_stat_sbn = cor(l2fc, sbn, method = "pearson"),
            cor_stat_spn = cor(l2fc, spn, method = "pearson"))
