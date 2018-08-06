library(tidyverse)


# Get families ------------------------------------------------------------

fams <- "https://funricegenes.github.io/famInfo.table.txt"

fams <- fams %>%
  read_delim(delim = "\t")



# From Nakano 2006 --------------------------------------------------------

ap2s <- "http://www.plantphysiol.org/highwire/filestream/120654/field_highwire_adjunct_files/1/73783supptable_III_rv_1202.xls"

dest <- "../data-raw/ap2s-nakano.xls"

if(!file.exists(dest)) ap2s %>% download.file(destfile = dest)

ap2s <- readxl::read_excel(path = dest, range = "A2:G141")

# Get ap2 subfams ---------------------------------------------------------

ap2_subfam <- c("PLT",
                "ERF",
                "DREB",
                "AP2")

fams %>%
  filter(Name %in% ap2_subfam)
