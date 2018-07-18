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


# Function that converts MSU to RAP ---------------------------------------

pull_rap <- function(ids) {
  id_rap <- dict %>%
    filter(grepl(paste(ids, collapse = "|"), msu)) %>%
    filter(rapdb != "None") %>%
    pull(rapdb)
  return(id_rap)
}

# Genes of SM and BM ------------------------------------------------------

load(file = "../data/rlog-pca.Rdata")

pc5 <- pcro %>%
  select(PC5, locus_id) %>%
  arrange(PC5)

last_pc5 <- pc5 %>%
  top_n(-100, wt = PC5) %>%
  pull(locus_id) 

top_pc5 <- pc5 %>%
  top_n(100, wt = PC5) %>%
  pull(locus_id)

last_pc5_400 <-  pc5 %>%
  top_n(-400, wt = PC5) %>%
  pull(locus_id)

top_pc5_400 <- pc5 %>%
  top_n(400, wt = PC5) %>%
  pull(locus_id)

# define function that downloads sequences --------------------------------

id <- "Os03g0215400"

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
                          '"format" : "fasta"}')) 
  if(status_code(r) != 200) {
    return(NA)
  } else {
    # ready for dataframe
    r <- r %>% httr::content() %>%
      purrr::flatten()
    
    r <- r[c("desc", "query", "id", "seq", "molecule")]
    return(r)
    }
}

# wrapper download sequences -------------------------------------------------------

ids <- top_pc5

get_top_seq <- function(ids) 
  {
  tst <- ids %>%
    pull_rap() %>%
    map(get_seq)
  
  sum(is.na(tst))
  
  n_pc5_df <- tst[!is.na(tst)] %>%
    purrr::reduce(.f = bind_rows) %>%
    mutate(id = paste(id, desc))
  
  path <- paste0("../seq/",
                 deparse(substitute(ids)),
                 "_2000TSS.fasta")
  
  set_names(x = n_pc5_df$seq,
            nm  = n_pc5_df$id) %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::subseq(start = 1, end = 2000) %>%
    writeXStringSet(filepath = path)
}


# Download sequences ------------------------------------------------------

get_top_seq(ids = top_pc5)

get_top_seq(ids = top_pc5_400)

get_top_seq(ids = last_pc5_400)

# get some random promoters as background
set.seed(1)
random_genes <- sample(pcro$locus_id, size = 400)

get_top_seq(random_genes)

a <- function(ids) {
  paste0("../seq/",
          deparse(substitute(ids)),
          "_2000TSS.fasta")
}

p <- "ciao"

b <- a(top_pc5)


# Get sequences of genes with nice patterns -------------------------------

nice_pattern <- c("Os10g0472900", "Os07g0136300", 
                  "Os04g0656500",
                  "Os10g0403000", "Os02g0492900", "Os03g0680200",
                  "Os01g0201600", "Os07g0683600", "Os01g0726400",
                  "Os04g0569100", "Os01g0832600")

# system("./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS.fasta \
#        -dna \
#        -nmotifs 10\
#        -p 4")

