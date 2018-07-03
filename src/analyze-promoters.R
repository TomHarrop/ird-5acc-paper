library(tidyverse)
library(httr)
library(jsonlite)
library(Biostrings)

# Download MSU to RAPDB ---------------------------------------------------
dict_path <- "../data/msu-to-rapdb.Rdata"


if(!file.exists(dict_path)) {
  dict_url <-  "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz"
  dict <- read_delim(dict_url,
                     delim = "\t", 
                     col_names = c("rapdb", "msu"))
  save(dict, file = dict_path)
} else {
  load(dict_path)
}  


# Functions that convert MSU tu RAP ---------------------------------------


id <- c("LOC_Os03g12950", "LOC_Os03g60430")
id_rap <- dict %>%
  filter(grepl(paste(id, collapse = "|"), msu)) %>%
  pull(rapdb)

pull_rap <- function(ids) {
  id_rap <- dict %>%
    filter(grepl(paste(ids, collapse = "|"), msu)) %>%
    pull(rapdb)
  return(id_rap)
}

# gramene / Ensembl REST server ------------------------------------------

server <- "http://rest.ensemblgenomes.org"
ext <- "/sequence/id"


# Define useful functions ----------------------------------------------------------

# this downloads promoter [2000 before TSS]
# and gene sequence

get_coords <- function(id) 
{
  # get coordinates
  r <- POST(paste(server, ext, sep = ""),
            content_type("application/json"),
            accept("application/json"),
            # body = '{ "ids" : ["Os03g0232200"],
            body = paste0('{ "ids" : ["', id,'"],',
                          '"expand_5prime" : [2000],',
          '"format" : "fasta"}')) %>%
    # At this point I just need the coordinates
    content() %>%
    purrr::flatten() %>%
    # purrr::reduce(., .f = bind_rows) %>%
    .$desc %>% 
    str_split_fixed(pattern = ":",
                    n = 3) %>% .[3] 
  
  # coords better stored in named character
  coords <- set_names(x = r, nm = id)
  
  return(coords)
}


# get synteny from coordinates in sativa

get_synteny <- function(region, species) {
  ext <- paste0("/alignment/region/oryza_sativa/",
                region, "?",
                "method=LASTZ_NET;",
                "aligned=0;",
                "species_set=", "oryza_sativa",";",
                # "species_set=oryza_indica;",
                "species_set=", species)  
  
  r <- GET(paste(server, ext, sep = ""),
           content_type("application/json"))
  
  return(r)
}



make_alignment_df <- function(alignments) {
  alig_df <- purrr:::reduce(alignments[[1]]$alignments,
                            .f = bind_rows) %>%
    mutate(names_desc = paste(species,seq_region,
                              start, end, strand, sep = ":"))
  return(alig_df)
}

# aligs <- map(coord, get_alignment) %>%
#   map(content)


# Wrapper that downloads synteny groups -----------------------------------

get_sequences <- function(coords, species) 
{
  # get syntenies for all species
  alig <- species %>%
    map(., ~get_synteny(region = coords,
                        species = .)) %>%
    map(content)
  
  # turn them into a tibble
  alig_df <- map(alig, make_alignment_df) %>%
    purrr::reduce(.f = bind_rows) %>%
    distinct()
  
  # and into a biostring object
  alig_df <- set_names(alig_df$seq, nm = alig_df$names_desc) %>%
    Biostrings::DNAStringSet()
}


# Test  workflow ----------------------------------------------------------

species <- c(barthii = "oryza_barthii",
             indica = "oryza_indica",
             rufipogon = "oryza_rufipogon",
             glaberrima = "oryza_glaberrima")


tst <- id_rap[2] %>%
  get_coords() %>%
  get_sequences(species = species) 
# %>%
#   subseq(1, 2000)

# Run on selected genes ---------------------------------------------------


"LOC_Os03g60430" %>%
  pull_rap() %>%
  get_coords() %>%
  get_sequences(species = species) %>%
  writeXStringSet(filepath = "../seq/LOC_Os03g60430.fasta")

# osZHD4 upregulated in SM in africans LOC_Os11g13930
# interesting, promoter not in synteny on african rice

"LOC_Os11g13930" %>%
  pull_rap() %>%
  get_coords() %>%
  get_sequences(species = species) %>%
  writeXStringSet(filepath = "../seq/LOC_Os11g13930.fasta")

# OsMADS4 - LOC_Os05g34940
# up in SM, expecially in african species (duh!)

osmads4 <- "LOC_Os05g34940" 

osmads4 %>%
  pull_rap() %>%
  get_coords() %>%
  get_sequences(species = species) %>%
  msa(method = "ClustalOmega") %>%
  DNAStringSet() %>%
  writeXStringSet(filepath = paste0("../seq/",
                                    osmads4, "-",
                                    quote(osmads4),
                                    ".fasta"))




# # with mysql
# 
# con <- dbConnect(RMySQL::MySQL(),
#                  host = "mysql-eg-publicsql.ebi.ac.uk",
#                  port = 4157,
#                  username = "anonymous")
# res <- dbSendQuery(con, "use family;")
# dbFetch(res)
# dbClearResult(res)
# 
# dbDisconnect(con)
# 
# 
# con <- dbConnect(RMySQL::MySQL(),
#                  host = "ensembldb.ensembl.org",
#                  port = 5306,
#                  username = "anonymous")
# dbGetInfo(con)
# rs <- dbSendQuery(con, "SHOW SCHEMAS LIKE '%compara%'")
# dbFetch(rs)
# dbClearResult(rs)
# dbDisconnect(con)
# con <- dbConnect(RMySQL::MySQL(),
#                  host = "ensembldb.ensembl.org",
#                  port = 5306,
#                  username = "anonymous", 
#                  dbname = "ensembl_compara_92")
# dbGetInfo(con)
# tst <- dbListTables(con)
# dbDisconnect(con)