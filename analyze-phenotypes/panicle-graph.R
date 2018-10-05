# install.packages("archiDART")
library(tidyverse)
library(magrittr)
library(xml2)
# library(tidygraph)
library(ggraph)
library(igraph)

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
  # geom_edge_link(arrow = grid::arrow(length = unit(0.08,
  #                                                  "inches"),
  #                                    type = "closed")) + 
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

# with Spatial and gDistance ----------------------------------------------

add_spikelet <- function(graph_in,
                         new_vert = numeric(2)) 
  {
  
  # new_vert in format c(x, y)
  
  # remove this stuff
  # new_vert_id <- 10
  # new_vert <- as.numeric(seeds[new_vert_id, 1:2])
  # graph_in <- tst
  
  # Long datagrame easier to subset than graph 
  graph_long <- graph_in %>% igraph::as_long_data_frame()
  
  # function
  get_distance <- function(x0, y0,
                           x1, y1,
                           x2, y2) 
    {
    ln <- sp::Line(matrix(c(x1, x2,
                            y1, y2),
                          ncol = 2)) %>%
      sp::Lines(., ID = "a") %>%
      list(.) %>%
      sp::SpatialLines()
    pn <- sp::SpatialPoints(coords = matrix(c(x0, y0),
                                            ncol = 2))
    # print(ln)
    # print(pn)
    return(rgeos::gDistance(ln, pn))
  }
  
  # test
  # get_distance(1328,  712, 1284,  688, 1329,  579)
  
  # Find edge rank
  nearest_edge <- graph_long %>%
    pmap(.,
         ~get_distance(x0 = new_vert[1],
                       y0 = new_vert[2],
                       x1 = ..3,
                       y1 = ..4,
                       x2 = ..8,
                       y2 = ..9)) %>%
    purrr::reduce(c) %>%
    which.min()
  
  # count edges
  vn <- vcount(graph_in)
  new_vn <- vn + 1
  
  # add spikelet and new edges
  new_graph <- graph_in %>%
    igraph::delete_edges(edges = paste0(as.numeric(graph_long[nearest_edge, "from"]),
                                        "|",
                                        as.numeric(graph_long[nearest_edge, "to"]))) %>%
    igraph::add_vertices(nv = 1,
                         attr = list(x = new_vert[1],
                                     y = new_vert[2],
                                     type = "spikelet")) %>%
    # igraph::as_long_data_frame()
    igraph::add_edges(edges = c(as.numeric(graph_long[nearest_edge, "from"]), new_vn,
                                new_vn, as.numeric(graph_long[nearest_edge, "to"])))
  
 
  return(new_graph) 
}

add_all_spikelets <- function(graph_base, spikelets)
{
  # spikelets <- seeds
  # graph_base <- tst
  assign_spikelet <- function(graph_base, spikelet) {
    graph_base <<- add_spikelet(graph_base, spikelet)
  }
  pmap(spikelets, ~assign_spikelet(graph_base = graph_base,
                                   spikelet = c(..1, ..2)))
  
  return(graph_base)
}

# Test function that adds spikelet ----------------------------------------

new_graph <- add_spikelet(graph_in = tst,
                          new_vert = as.numeric(seeds[10, ]))

new_graph %>% ggraph() +
  geom_edge_link() +
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()
  
all_spikelets <- add_all_spikelets(graph_base = tst,
                                   spikelets = seeds)

all_spikelets %>%
  ggraph() +
  # geom_edge_link(arrow = grid::arrow(length = unit(0.08,
  #                                                  "inches"),
  #                                    type = "closed")) +
  geom_edge_link() +
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()

