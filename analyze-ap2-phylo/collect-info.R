library(tidyverse)
library(magrittr)

# Get families ------------------------------------------------------------

fams <- "https://funricegenes.github.io/famInfo.table.txt"

fams <- fams %>%
  read_delim(delim = "\t") %>%
  mutate(locus_id = MSU)



# From Nakano 2006 --------------------------------------------------------

ap2s <- "http://www.plantphysiol.org/highwire/filestream/120654/field_highwire_adjunct_files/1/73783supptable_III_rv_1202.xls"

dest <- "../data-raw/ap2s-nakano.xls"

if(!file.exists(dest)) ap2s %>% download.file(destfile = dest)

ap2s <- readxl::read_excel(path = dest, range = "A2:G141") %>%
  mutate(locus_id = paste0("LOC_", `TIGR Locus identifier`))

# Get ap2 subfams ---------------------------------------------------------

ap2_subfam <- c("PLT",
                "ERF",
                "DREB",
                "AP2")

fams_erf <- fams %>%
  filter(Name %in% ap2_subfam)

all_ap2s <- fams_erf %>% 
  full_join(ap2s) %>%
  write_csv("../data/ap2_full.csv")


# Plot all ap2s -----------------------------------------------------------

all_ap2s <- all_ap2s %>% 
  select(locus_id, Name, `Generic Name`) %>%
  dplyr::rename(generic_name = `Generic Name`)

library(DESeq2)

source("../src/helper-functions.R")
dds <- readRDS(file = "../data-raw/dds.Rds")

load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
load("../data/msu-to-rapdb.Rdata"); rap2msu <- dict; rm(dict)

tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

annos <- annos %>%
  dplyr::rename(locus_id = "MSU") %>%
  select(symbol, locus_id)

mapman <- mapman %>%
  dplyr::rename(locus_id = "IDENTIFIER") %>%
  filter(!duplicated(locus_id))



pdf("../fig/more_ap2s.pdf", height = 80, width = 20)
all_ap2s %$%
  get_expression(locus_id, dds = dds) %>%
  left_join(all_ap2s) %>%
  left_join(annos) %>%
  left_join(tf_fam) %>%
  left_join(mapman) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id",
                        "Name",
                        "generic_name",
                        "symbol",
                        "Family",
                        "DESCRIPTION"),
             scales = "free_y",
             ncol = 5,
             labeller = label_wrap_gen(width = 50,
                                       multi_line = T))
dev.off()  
