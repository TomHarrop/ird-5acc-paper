# install.packages("archiDART")
library(tidyverse)
library(magrittr)
library(xml2)
# library(tidygraph)
library(ggraph)

vertices <- read_xml("../data-raw/pheno-xml/Nip_1_1_6307.ricepr") %>%
  # xml_child("graph")  %>%
  # xml_child("vertices") %>%
  xml_find_all(".//vertex")
  
nodes <- tibble(x = as.numeric(xml_attr(vertices, "x")),
       y = as.numeric(xml_attr(vertices, "y")),
       type = xml_attr(vertices, "type"),
       id = xml_attr(vertices, "id"))

edges <- read_xml("../data-raw/pheno-xml/Nip_1_1_6307.ricepr") %>%
  xml_find_all(".//edge")

edges <- tibble(from = xml_attr(edges, "vertex1"),
                to = xml_attr(edges, "vertex2"))

tst <- tidygraph::tbl_graph(nodes = nodes,
                            edges = edges, )
  
ggraph(tst)+ 
  geom_edge_link() + 
  geom_node_point(size = 8, colour = 'steelblue') 

ggplot(dat,
       aes(x = x,
           y = y,
           colour = type)) +
  geom_point() +
  theme_bw()

# check
# https://www.data-imaginist.com/2017/introducing-tidygraph/
