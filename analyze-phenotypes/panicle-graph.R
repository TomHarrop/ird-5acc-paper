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


# Spatial with gDistance --------------------------------------------------

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
  print(ln)
  print(pn)
  rgeos::gDistance(ln, pn)
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


# or ----------------------------------------------------------------------

# from https://stackoverflow.com/questions/35194048/
# using-r-how-to-calculate-the-distance-from-one-point-to-a-line
# https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
# dist2d <- function(a,b,c) {
#   print(c(a, b, c))
#   v1 <- b - c
#   v2 <- a - b
#   m <- cbind(v1,v2)
#   d <- abs(det(m))/sqrt(sum(v1*v1))
# } 

# recode function
eu_dist <- function(p0, p1, p2) {
  print(c(p0, p1, p2))
  # p0 is the point
  # the line goes through p1 and p1
  # every point is p <- c(x, y)
  N <- abs(((p2[2] - p1[2])*p0[1]) -
             ((p2[1] - p1[1])*p0[2]) +
             (p2[1]*p1[2]) -
             (p2[2]*p1[1]))
  D <- sqrt(((p2[2] - p1[2])^2) +
              ((p2[1] - p1[1])^2))
  print(N)
  print(D)
  distance <- N/D
  print(distance)
  return(N/D)
}


tst_df <- tst %>% igraph::as_long_data_frame()

new_vert <- 10
new_vert_coord <- as.numeric(seeds[new_vert, 1:2])

# nearest_edge <- tst_df %>%
#   pmap(.,
#        ~dist2d(a = as.numeric(seeds[new_vert, 1:2]),
#                b = c(..3, ..4),
#                c = c(..8, ..9))) %>%
#   purrr::reduce(c) %>%
#   which.min()

nearest_edge <- tst_df %>%
  pmap(.,
       ~eu_dist(p0 = as.numeric(seeds[new_vert, 1:2]),
                p1 = c(..3, ..4),
                p2 = c(..8, ..9))) %>%
  purrr::reduce(c) %>%
  which.min()


# this is not vectorized
# nearest_edge <- tst_df %>% 
#   mutate(dist = dist2d(a = c(1328,  712),
#                        b = c(from_x, from_y),
#                        c = c(to_x, to_y))) %>%
#   purrr::reduce(c) %>%
#   which.min()


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

