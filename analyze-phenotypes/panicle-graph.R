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
       id = xml_attr(vertices, "id")) %>%
  mutate(rank = 1:n())

nodes_rank <- set_names(x = nodes$rank, nm = nodes$id)

edges <- read_xml("../data-raw/pheno-xml/Nip_1_1_6307.ricepr") %>%
  xml_find_all(".//edge")

edges <- tibble(from = xml_attr(edges, "vertex1"),
                to = xml_attr(edges, "vertex2")) %>%
  mutate(from = nodes_rank[from],
         to = nodes_rank[to]) %>%
  mutate_all(unname)

tst <- tidygraph::tbl_graph(nodes = nodes,
                            edges = edges)
  
ggraph(tst) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()

# ggplot(dat,
#        aes(x = x,
#            y = y,
#            colour = type)) +
#   geom_point() +
#   theme_bw()


# Read XML file with seeds ------------------------------------------------

seeds <- read_xml("../data-raw/pheno-xml/Nip_1_1_6307.ricegr") %>%
  xml_find_all(".//particle") 

seeds <- tibble(x = as.numeric(xml_attr(seeds, "cx")),
                y = as.numeric(xml_attr(seeds, "cy")))

ggraph(tst) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type),
                  size = 2) +
  geom_point(data = seeds, 
             aes(x = x, 
                 y = y)) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()
# check
# https://www.data-imaginist.com/2017/introducing-tidygraph/

# how to guess path
# spatial network
# spatial network add vertex to edge
# points2network() ?