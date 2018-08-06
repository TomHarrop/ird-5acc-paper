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

pdf("../fig/more_ap2s.pdf", height = 50, width = 20)
all_ap2s %$%
  get_expression(locus_id, dds = dds) %>%
  left_join(all_ap2s) %>%
  plot_norm_expr() +
  facet_wrap(facets = c("locus_id", "Name", "generic_name"),
             scales = "free_y",
             ncol = 5)
dev.off()  
