library(tidyverse)

load("../data/all-sep-deseq.Rdata")

spc_res <- spc_res %>% map(~as.data.frame(.[, c("stat", "locus_id")]))

rename_stat <- function(i) {
  colnames(spc_res[[i]])[1] <- paste("stat", i, sep = "_")
  return(spc_res[[i]])
}

tst <- names(spc_res) %>% map(rename_stat)

tst <- tst %>% purrr::reduce(full_join, by = "locus_id")

rownames(tst) <- tst$locus_id

tst$locus_id <- NULL
tst[is.na(tst)] <- 0

pc <- prcomp(tst, scale. = T)

pc$sdev

biplot(pc)
