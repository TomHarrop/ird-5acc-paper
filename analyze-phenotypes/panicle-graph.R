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


# Test Spatial --------------------------------------------------------------
# SpatialLinesDataFrame 

# tst_df <- tst %>%
#   igraph::as_long_data_frame() %>%
#   remove_rownames() %>%
#   column_to_rownames("to_rank")
# 
# tst_sp <-  tst %>%
#   igraph::as_long_data_frame() %>%
#   pmap(., ~sp::Line(matrix(c(..3, ..8, ..4, ..9), ncol = 2)) %>%
#          sp::Lines(ID = ..12)) %>%
#   sp::SpatialLines() %>%
#   sp::SpatialLinesDataFrame(data = tst_df, ) %T>%
#   plot()
# 
# library(shp2graph)
# tst_seeds <- points2network(ntdata = tst_sp,
#                           pointsxy = seeds) %>%
#   shp2graph::nel2igraph()
# 
# # from http://r-sig-geo.2731867.n2.nabble.com/igraph-and-spatial-td7589564.html
# makeLineFromCoords <- function(coords, i) { 
#   Sl1 = Line(coords) 
#   S1 = Lines(list(Sl1), ID=as.character(i)) 
#   Sl = SpatialLines(list(S1)) 
#   return(Sl) 
# } 

# overcomplicated for what we have to achieve

# or ----------------------------------------------------------------------

# from https://stackoverflow.com/questions/35194048/
# using-r-how-to-calculate-the-distance-from-one-point-to-a-line
# https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
dist2d <- function(a,b,c) {
  # print(c(a, b, c))
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 

tst_df <- tst %>% igraph::as_long_data_frame()

new_vert <- 10

tst_df %>% pmap(., 
                ~dist2d(a = as.numeric(seeds[new_vert, 1:2]),
                        b = c(..3, ..4),
                        c = c(..8, ..9))) %>%
  purrr::reduce(c) %>%
  which.min()


new_graph <- tst %>% 
  igraph::add_vertices(nv = 1, 
                       attr = list(x = as.numeric(seeds[new_vert, 1]),
                                   y = as.numeric(seeds[new_vert, 2]),
                                   type = "spikelet")) %>%
  # igraph::as_long_data_frame()
  igraph::add_edges(edges = c(as.numeric(tst_df[20, "from"]), 50,
                              50, as.numeric(tst_df[20, "to"])))


ggraph(new_graph) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()

