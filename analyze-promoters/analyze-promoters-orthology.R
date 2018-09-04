library(tidyverse)
library(httr)
library(jsonlite)
library(Biostrings)
library(msa)

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


# define species ----------------------------------------------------------

ory_species <- c(barthii = "oryza_barthii",
                 indica = "oryza_indica",
                 rufipogon = "oryza_rufipogon",
                 glaberrima = "oryza_glaberrima")


# function that fetch orthologs ---------------------------------------------------------

fetch_ortho <- function(id) {
  server <- "http://rest.ensemblgenomes.org"
  
  ext <- paste0("/homology/id/",
                id,
                "?")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)

  
  r <- head(fromJSON(toJSON(content(r)))) %>%
    purrr::flatten() %>%
    purrr::flatten() %>%
    purrr::flatten() 
    
  r <- r$target %>%
    map(unlist) %>%
    do.call(what = bind_cols) %>%
    filter(species %in% ory_species)
  
  return(r)  
}

# Define function that fetch sequence -------------------------------------

id <- "Os03g0232200"
get_seq <- function(id) 
{
  server <- "http://rest.ensemblgenomes.org"
  ext <- "/sequence/id"
  
  # get coordinates
  r <- POST(paste(server, ext, sep = ""),
            content_type("application/json"),
            accept("application/json"),
            # body = '{ "ids" : ["Os03g0232200"],
            body = paste0('{ "ids" : ["', id,'"],',
                          '"expand_5prime" : [2000],',
                          '"format" : "fasta"}')) %>%
    # ready for dataframe
    content() %>%
    purrr::flatten()
    
  return(r)
}

get_seq("Os03g0232200")


# Wrapper for the full function -------------------------------------------

# some_ap2 <- "LOC_Os03g12950" %>% pull_rap()

get_orthoseq <- function(id) 
{
  # get sequences
  r <- id %>%
    # get ids
    fetch_ortho() %>%
    .$id %>%
    c(id) %>%
    map(get_seq) 
  
  # turn seqs into data frame
  r <- r %>%
    map(purrr::flatten) %>%
    purrr::reduce(bind_rows) %>%
    mutate(desc = paste(id, desc))
  
  # make a biostring element
  r <- set_names(x = r$seq, r$desc) %>%
    Biostrings::DNAStringSet() %>%
    msa(method = "ClustalOmega") %>%
    DNAStringSet() 
  
  return(r)
}


# check promoters ---------------------------------------------------------

make_path <- function(gene) {
  paste0("../seq/",
         gene, "-",
         substitute(gene), "-",
         "al",
         ".fasta")
}


some_ap2 <- "LOC_Os03g12950"

some_ap2 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(some_ap2))

rav2 <- "LOC_Os01g04800"

rav2 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(rav2))


# OsMADS4 - LOC_Os05g34940
# up in SM, expecially in african species (duh!)

osmads4 <- "LOC_Os05g34940" 

osmads4 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(osmads4))

# osZHD4 upregulated in SM in africans LOC_Os11g13930
# interesting, promoter not in synteny on african rice

osZHD4 <- "LOC_Os11g13930"

osZHD4 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(osZHD4))


# osmads16-spw1 upregulated in BM

spw1 <- "LOC_Os06g49840"

spw1 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(spw1))

# Homeobox LOC_Os10g26500 
# up in BM in african species

hb_something <- "LOC_Os10g26500"

hb_something %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(hb_something))

# LSY1 LOC_Os05g02730
# up in SM

lsy1 <- "LOC_Os05g02730"

lsy1 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(lsy1))

# LOC_Os06g44860
# SPL10 according to some annotation
# up in the sm

maybe_spl10 <- "LOC_Os06g44860"

maybe_spl10 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(maybe_spl10))

# Growth factors
# cluster 2 domestication

# LOC_Os03g47140	OsGRF9

os_grf9 <- "LOC_Os03g47140"

os_grf9 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(os_grf9))

# LOC_Os03g51970	OsGRF6

os_grf6 <- "LOC_Os03g51970"

os_grf6 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(os_grf6))


# LOC_Os06g02560	OsGRF5

os_grf5 <- "LOC_Os06g02560"

os_grf5 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(os_grf5))


# GRF in Cluster 2 --------------------------------------------------------

cluster2_grf9 <- "LOC_Os03g47140"

cluster2_grf9 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(cluster2_grf9))

cluster2_grf6 <- "LOC_Os03g51970"

cluster2_grf6 %>%
  pull_rap() %>%
  get_orthoseq() %>%
  writeXStringSet(filepath = make_path(cluster2_grf6))


# alignments for helen ----------------------------------------------------

make_path2 <- function(gene) {
  paste0("../seq/",
         gene, "-",
         "al",
         ".fasta")
}

to_align_df <- readxl::read_excel("../data-raw/helen-alignment.xlsx") %>% 
  mutate(object = LOC_Name,
         nm = paste(LOC_Name, X__1, "cluster", Cluster, sep = "_")) %>%
  mutate(nm = str_replace(nm, "/", "_")) %>%
  mutate(nm = str_replace(nm, " ", "_"))
 
to_align <- setNames(object = to_align_df$object,
           nm = to_align_df$nm)

get_al_seq <- function(nm) {
  print(nm)
  pull_rap(to_align[[nm]]) %>%
    get_orthoseq() %>%
    writeXStringSet(filepath = make_path2(as.character(nm)))
}

names(to_align) %>% 
  map(get_al_seq)

# ISSUE with OsGR5
# "LOC_Os06g02560_OsGRF5_cluster_2" %>% map(get_al_seq)
# to_align["LOC_Os06g02560_OsGRF5_cluster_2"]

which(names(to_align) == "LOC_Os06g02560_OsGRF5_cluster_2")

to_align <- to_align[(which(names(to_align) == "LOC_Os06g02560_OsGRF5_cluster_2") + 1):length(to_align)]

names(to_align) %>% 
  map(get_al_seq)
# fetch orthologs method 2 ------------------------------------------------

# 
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# server <- "http://rest.ensemblgenomes.org"
# ext <- "/homology/id/AT3G52430?compara=plants"
# 
# r <- GET(paste(server, ext, sep = ""), content_type("text/xml"))
# 
# stop_for_status(r)
# 
# 
# print(content(r))
/