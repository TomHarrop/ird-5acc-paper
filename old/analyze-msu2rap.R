library(tidyverse)

load("../data/msu-to-rapdb.Rdata")

dds <- readRDS(file = "../data-raw/dds.Rds")

locus_id <- rownames(dds)

# locus_id <- sample(locus_id, 100)

locus_id <- set_names(x = locus_id, nm = locus_id)

get_rap_ids <- function(id) 
{
  rap_id <- dict %>%
    filter(grepl(id, msu)) %>%
    pull(rapdb)
  
  return(paste(rap_id, collapse = ", "))
}

msu2rap <- map(locus_id, get_rap_ids)
msu2rap <- tibble(locus_id = names(msu2rap),
                  rap_id = flatten_chr(msu2rap))

save(msu2rap, file = "../data/msu2rap.Rdata")
