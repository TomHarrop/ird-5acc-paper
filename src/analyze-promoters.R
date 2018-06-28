library(tidyverse)
library(httr)
library(jsonlite)
library(RMySQL)

# Download MSU to RAPDB ---------------------------------------------------

dict_url <-  "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz"


dict_path <- "../data/msu-to-rapdb.Rdata"
if(!file.exists(dict_path)) {
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

ext <- paste0("/sequence/id/", "LOC_Os12g36040.1",
              "?object_type=transcript;",
              "db_type=otherfeatures;",
              "type=cds;",
              "species=oryza_sativa")



server <- "http://rest.ensemblgenomes.org"
ext <- "/sequence/id"
r <- POST(paste(server, ext, sep = ""),
          content_type("application/json"),
          accept("application/json"),
          # body = '{ "ids" : ["Os03g0232200", "Os03g0818800"],
          body = '{ "ids" : ["Os03g0232200"],
          "format" : "fasta"}')

content(r)

# this downloads promoter [2000 before TSS]
# and gene sequence

r_ext <- POST(paste(server, ext, sep = ""),
              content_type("application/json"),
              accept("application/json"),
              # body = '{ "ids" : ["Os03g0232200"],
              body = '{ "ids" : ["Os03g0232200", "Os03g0818800"],
                        "expand_5prime" : [2000],
          "format" : "fasta"}')

tst <- content(r_ext) %>%
  purrr::reduce(., .f = bind_rows) %>%
  mutate(coord = str_split_fixed(string = desc,
                                 pattern = ":",
                                 n = 3)[,3])
  

# write(content(r), file = "tst.fasta")

Biostrings::DNAStringSet(x = tst$seq)



# Get alignments ----------------------------------------------------------

server <- "http://rest.ensemblgenomes.org"

get_alignment <- function(region) {
  ext <- paste0("/alignment/region/oryza_sativa/",
                region, "?",
                "method=LASTZ_NET;",
                "species_set=", "oryza_sativa",";",
                # "species_set=oryza_indica;",
                "species_set=oryza_barthii")  
  
  r <- GET(paste(server, ext, sep = ""),
           content_type("application/json"))
  
  return(r)
}

tst$coord %>% map(get_alignment)

r <- GET(paste(server, ext, sep = ""),
         content_type("application/json"))


stop_for_status(r)
content(r)

# with mysql

con <- dbConnect(RMySQL::MySQL(),
                 host = "mysql-eg-publicsql.ebi.ac.uk",
                 port = 4157,
                 username = "anonymous")
res <- dbSendQuery(con, "use family;")
dbFetch(res)
dbClearResult(res)

dbDisconnect(con)


con <- dbConnect(RMySQL::MySQL(),
                 host = "ensembldb.ensembl.org",
                 port = 5306,
                 username = "anonymous")
dbGetInfo(con)
rs <- dbSendQuery(con, "SHOW SCHEMAS LIKE '%compara%'")
dbFetch(rs)
dbClearResult(rs)
dbDisconnect(con)
con <- dbConnect(RMySQL::MySQL(),
                 host = "ensembldb.ensembl.org",
                 port = 5306,
                 username = "anonymous", 
                 dbname = "ensembl_compara_92")
dbGetInfo(con)
tst <- dbListTables(con)
dbDisconnect(con)