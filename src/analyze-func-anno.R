library(tidyverse)


# if(!file.exists("../data/func-anno.Rdata")) {
#   urls <- list(gene_info = "https://funricegenes.github.io/geneInfo.table.txt",
#                gene_family = "https://funricegenes.github.io/famInfo.table.txt",
#                # papers = "https://funricegenes.github.io/reference.table.txt",
#                keywords = "https://funricegenes.github.io/geneKeyword.table.txt")
#   
#   annos <- urls %>% map(~read_delim(., delim = "\t"))
#   save(annos, file = "../data/func-anno.Rdata")
# }
# 
# load("../data/func-anno.Rdata")

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

# library(tidyverse)
# 
# rapd_msu_file <- "../data-raw/RAP-MSU_2018-03-29.txt.gz" 
# rapd_anno_file <- "../data-raw/IRGSP-1.0_predicted_annotation_2018-03-29.tsv.gz"
# 
# 
# if(!file.exists(rapd_msu_file)) {
# download.file("http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz",
#               destfile = rapd_msu_file)
# }
# 
# if(!file.exists(rapd_anno_file)) {
#   download.file("http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_predicted_annotation_2018-03-29.tsv.gz",
#                 destfile = rapd_anno_file)
# }
# 
# rapd_msu <- read_delim(rapd_msu_file,
#                        delim = "\t", col_names = c("rapdb", "msu")) 
# 
# 
# rapd_anno <- read_delim(rapd_anno_file, delim = "\t") %>%
#   select(Locus_ID, `CGSNL Gene Symbol`, `CGSNL Gene Name`, `Oryzabase Gene Symbol Synonym(s)`) %>%
#   rename(Locus_ID = "rapdb")
# 
# 
# # There is more than one msu id for RAP id (also after removing splice variants)
# # tst_ <- rapd_msu %>%
# #   mutate(msu_s = str_split(msu, pattern = ",")) %>%
# #   mutate(msu_s = map(msu_s, ~ str_split(., pattern = "\\.", simplify = TRUE))) %>%
# #   mutate(len = map(msu_s, ~length(unique(.[, 1])))) %>%
# #   mutate(len = unlist(len)) %>%
# #   filter(len > 1)
# # 
# # summary(tst$len)
# 
# # How to deal with this?
# 
# rapd_msu <- rapd_msu %>%
#   mutate(msu = str_split(msu, pattern = ",")) %>%
#   mutate(msu = map(msu, ~ str_split(., pattern = "\\.", simplify = TRUE)[, 1])) %>%
#   mutate(msu = map(msu, unique)) %>% 
#   unnest() %>%
#   distinct()
# 
# # match msu id with annotations
# 
# tst <- rapd_msu %>% 
#   full_join(rapd_anno) %>%
#   filter(!is.na(`Oryzabase Gene Symbol Synonym(s)`))



