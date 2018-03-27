# Set up environment ------------------------------------------------------

library(tidyverse)
library(readxl)

# Load Mapman
# 
# requires readxl and tidyverse package
# 
# The Load Mapman function load and wrangles mapman files so that they can
# be used for testing category enrichement
# 
# @param file a mapman file in the excel format
# 
# Test
# mapman <- load_mapman(file = "data/MAPMAN BIN-Osa_MSU_v7.xlsx")

load_mapman <- function(path)
{
  mapman <- read_excel(path)
  mapman <- mapman %>% filter(! is.na(IDENTIFIER))
  
  # How many levels in the bins?
  max_l <- max(sapply(str_split(mapman$NAME, "\\."), length))
  
  # So, make it a data frame
  mapman_names <- str_split_fixed(mapman$NAME, "\\.", max_l) 
  mapman_names[nchar(mapman_names) == 0] <- NA 
  rownames(mapman_names) <- mapman$IDENTIFIER
  
  # you need some sort of cumulative string split
  # This is a very inelegant solution
  mapman_bins <- seq_along(rownames(mapman_names)) %>% #t(mapman_names)
    map(~map_chr(1:7, function(n) {
      gsub(".NA", "",
           paste(mapman_names[., 1:n], collapse = "."),
           fixed = TRUE)
    })) %>%
    do.call(rbind, .)
  mapman_bins[is.na(mapman_names)] <- NA
  mapman_bins <- data.frame(mapman_bins, stringsAsFactors = FALSE)
  colnames(mapman_bins) <- paste0("level", 1:ncol(mapman_bins))
  mapman_bins$IDENTIFIER <- mapman$IDENTIFIER
  
  # merge back
  # I would use full_join(), but IDENTIFIER are dublicated
  mapman <- bind_cols(mapman, mapman_bins)
  
  return(mapman)
}



mapman <- load_mapman("../data/MAPMAN BIN-Osa_MSU_v7.xlsx")
save(mapman, file = "../data/mapman.Rdata")
