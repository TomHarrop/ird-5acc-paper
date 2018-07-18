library(tidyverse)

paths <- list(gene_info = "../data-raw/geneInfo.table.txt",
             gene_family = "../data-raw/famInfo.table.txt",
             # papers = "https://funricegenes.github.io/reference.table.txt",
             keywords = "../data-raw/geneKeyword.table.txt")

annos <- paths %>% map(~read_delim(., delim = "\t"))
# save(annos, file = "../data/func-anno.Rdata")
  
  
annos <- full_join(annos$gene_info, annos$gene_family, by = 'MSU') %>%
  filter(MSU != "None") %>%
  rename(symbol = "Symbol.x",
         symbol_from_family = "Symbol.y",
         name_of_family  = "Name") %>%
  full_join(annos$keywords) %>%
  distinct() %>%
  filter(MSU != "None") %>%
  group_by_at(vars(symbol:RAPdb)) %>%
  # nest() %>%
  summarise(Keyword = paste(Keyword, collapse = " - "),
            Title = paste(Title, collapse = " --- ")) %>%
  ungroup() %>%
  # mutate(RAPdb = list(RAPdb.x, RAPdb.y, RAPdb))
  mutate(RAPdb = paste(RAPdb.x, RAPdb.y, RAPdb)) %>%
  mutate(RAPdb = map(RAPdb, ~str_split(., pattern = " ", simplify = T)[[1]][1])) %>%
  mutate(RAPdb = unlist(RAPdb)) %>%
  # mutate(RAPdb = map(RAPdb, ~unique(.[.!="NA"]))) %>%
  select(-RAPdb.x, -RAPdb.y) %>%
  rename(symbol_from_paper = "Symbol")
  
save(annos, file = "../data/func-anno.Rdata")
write.csv2(annos, "../data/annos_funricegene.csv") 


