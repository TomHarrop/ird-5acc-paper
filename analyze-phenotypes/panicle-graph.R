# install.packages("archiDART")
library(tidyverse)
library(magrittr)
library(xml2)

vertices <- read_xml("../data-raw/pheno-xml/Nip_1_1_6307.ricepr") %>%
  # xml_child("graph")  %>%
  # xml_child("vertices") %>%
  xml_find_all(".//vertex")
  
dat <- tibble(x = as.numeric(xml_attr(vertices, "x")),
       y = as.numeric(xml_attr(vertices, "y")),
       type = xml_attr(vertices, "type"))

ggplot(dat,
       aes(x = x,
           y = y,
           colour = type)) +
  geom_point() +
  theme_bw()

# check
# https://www.data-imaginist.com/2017/introducing-tidygraph/