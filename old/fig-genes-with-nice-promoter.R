library(tidyverse)

# Load datasets -----------------------------------------------------------

source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

load("../data/rlog-pca.Rdata")
load("../data/func-anno.Rdata")
load("../data/mapman.Rdata")
load("../data/msu2rap.Rdata")


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


# Get pc5 -----------------------------------------------------------------

pc5 <- pcro %>%
  select(PC5, locus_id)

last_pc5 <- pc5 %>%
  top_n(-100, wt = PC5) %>%
  pull(locus_id) 

top_pc5 <- pc5 %>%
  top_n(100, wt = PC5) %>%
  pull(locus_id)


# Plot --------------------------------------------------------------------
plot_ids <- function(ids, file_path) 
{
  p <- ids %>%
    get_expression(dds = dds) %>%
    left_join(annos) %>%
    left_join(tf_fam) %>%
    left_join(mapman) %>%
    left_join(msu2rap) %>%
    mutate(locus_id = as_factor(locus_id)) %>%
    plot_norm_expr() +
    facet_wrap(facets = c("locus_id",
                          "symbol",
                          "Family",
                          "DESCRIPTION",
                          "rap_id"),
               scales = "free_y",
               ncol = 5,
               labeller = label_wrap_gen(width = 50,
                                         multi_line = T))
  
  nplots <- length(levels(p$data$locus_id))
  
  pdf(file = file_path,
      height = ceiling(nplots/5)*5,
      width = ifelse(nplots < 5, nplots * 4, 20))
  print(p)
  dev.off()
}

nice_pattern <- c("Os10g0472900", "Os07g0136300", 
                  "Os04g0656500",
                  "Os10g0403000", "Os02g0492900", "Os03g0680200",
                  "Os01g0201600", "Os07g0683600", "Os01g0726400",
                  "Os04g0569100", "Os01g0832600")

msu2rap %>%
  filter(grepl(paste(nice_pattern, collapse = "|"),
               rap_id)) %>%
  pull(locus_id) %>%
  plot_ids("../fig/fig-genes-with-nice-promoter.pdf")



# Read MAST output --------------------------------------------------------

library(xml2)

ids_top <- xml2::read_xml("../seq/top-pc5-mast-out/mast.xml") %>%
  xml_child("sequences") %>%
  xml_contents() %>%
  xml_attr(attr = "name")

which(ids_top == "Os08g0159500")
nice_pattern_top400 <- ids_top[1:43]

msu2rap %>%
  filter(grepl(paste(nice_pattern_top400, collapse = "|"),
               rap_id)) %>%
  pull(locus_id) %>%
  plot_ids(file_path = "../fig/fig-genes-with-nice-promoter_top_400.pdf")

#########

ids_last <- xml2::read_xml("../seq/last-pc5-mast-out/mast.xml") %>%
  xml_child("sequences") %>%
  xml_contents() %>%
  xml_attr(attr = "name")

which(ids_last == "Os05g0513100")
nice_pattern_last400 <- ids_last[1:53]

msu2rap %>%
  filter(grepl(paste(nice_pattern_last400, collapse = "|"),
               rap_id)) %>%
  pull(locus_id) %>%
  plot_ids(file_path = "../fig/fig-genes-with-nice-promoter_last_400.pdf")


