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

tst_df <- tst %>% igraph::as_long_data_frame()
new_vert <- 10
new_vert_coord <- as.numeric(seeds[new_vert, 1:2])

get_distance <- function(x0, y0,
                         x1, y1,
                         x2, y2) {
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

get_distance(1328,  712, 1284,  688, 1329,  579)

nearest_edge <- tst_df %>%
  pmap(.,
       ~get_distance(x0 = as.numeric(seeds[new_vert, 1]),
                     y0 = as.numeric(seeds[new_vert, 2]),
                     x1 = ..3,
                     y1 = ..4,
                     x2 = ..8,
                     y2 = ..9)) %>%
  purrr::reduce(c) %>%
  which.min()

edge_out <- tst %>%
  igraph::add_vertices(nv = 1, 
                       attr = list(x = as.numeric(seeds[new_vert, 1]),
                                   y = as.numeric(seeds[new_vert, 2]),
                                   type = "spikelet")) %>%
  igraph::delete_edges(nearest_edge)

edge_out %>%
  ggraph() + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()

# count edges
vn <- vcount(tst)
new_vn <- vn + 1

new_graph <- tst %>%
  igraph::delete_edges(edges = paste0(as.numeric(tst_df[nearest_edge, "from"]),
                                      "|",
                                      as.numeric(tst_df[nearest_edge, "to"]))) %>%
  igraph::add_vertices(nv = 1,
                       attr = list(x = as.numeric(seeds[new_vert, 1]),
                                   y = as.numeric(seeds[new_vert, 2]),
                                   type = "spikelet")) %>%
  # igraph::as_long_data_frame()
  igraph::add_edges(edges = c(as.numeric(tst_df[nearest_edge, "from"]), new_vn,
                              new_vn, as.numeric(tst_df[nearest_edge, "to"])))

new_graph %>% ggraph() +
  geom_edge_link() +
  geom_node_point(aes(colour = type),
                  size = 2) +
  coord_fixed() +
  # coord_flip() +
  theme_minimal()
  
