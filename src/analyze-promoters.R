library(tidyverse)
library(httr)
library(jsonlite)
# library(RMySQL)

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

id <- c("LOC_Os03g12950", "LOC_Os03g60430")
id_rap <- dict %>%
  filter(grepl(paste(id, collapse = "|"), msu)) %>%
  pull(rapdb)


# connect to gramene REST server ------------------------------------------

server <- "http://rest.ensemblgenomes.org"
ext <- "/sequence/id"

# this downloads promoter [2000 before TSS]
# and gene sequence
r <- POST(paste(server, ext, sep = ""),
              content_type("application/json"),
              accept("application/json"),
              # body = '{ "ids" : ["Os03g0232200"],
              body = '{ "ids" : ["Os03g0232200", "Os03g0818800"],
                        "expand_5prime" : [2000],
          "format" : "fasta"}') %>%
  # At this point I just need the coordinates
  content() %>%
  purrr::reduce(., .f = bind_rows) %>%
  mutate(coord = str_split_fixed(string = desc,
                                 pattern = ":",
                                 n = 3)[,3]) 

coord <- set_names(x = r$coord, nm = r$id)

# Get alignments ----------------------------------------------------------

# get synteny from coordinates in sativa
get_alignment <- function(region, species) {
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

# aligs <- map(coord, get_alignment) %>%
#   map(content)

species <- c(barthii = "oryza_barthii",
             indica = "oryza_indica",
             rufipogon = "oryza_rufipogon",
             glaberrima = "oryza_glaberrima")

alig <- species %>%
  map(., ~get_alignment(region = coord[2],
                     species = .)) %>%
  map(content)

make_alignment_df <- function(alignments) {
  alig_df <- purrr:::reduce(alignments[[1]]$alignments,
                            .f = bind_rows) %>%
    mutate(names_desc = paste(species,seq_region,
                              start, end, strand, sep = ":"))
  return(alig_df)
}

alig_df <- map(alig, make_alignment_df) %>%
  purrr::reduce(.f = bind_rows) %>%
  distinct()

alig_df <- set_names(alig_df$seq, nm = alig_df$names_desc) %>%
  Biostrings::DNAStringSet()

Biostrings::writeXStringSet(alig_df, filepath = "../seq/tst.fasta")

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