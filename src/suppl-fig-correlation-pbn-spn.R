library(tidyverse)
library(ggfortify)
library(gridExtra)

load("../data/phenotypes.Rdata")

# Function that plots correlation -----------------------------------------


plot_pheno <- function(data = pheno_cali,
                       x = pbn,
                       y = spn)
{
  x <- enquo(x)
  y <- enquo(y)
  
  p <- data %>% 
    # select(!!x)
    group_by(Species) %>%
    mutate(cors = cor(!!x, !!y) %>%
             round(2) %>%
             {paste0("r = ", .)}) %>%
    ungroup()  %>%
    mutate(Species = paste(Species, cors, sep = " | ")) %>%
    ggplot(aes(x = !!x, y = !!y)) +
    # geom_point(alpha = .2) +
    geom_hex(bins = 20) +
    facet_wrap(facets = "Species",
               nrow = 1) +
    theme_bw() +
    scale_fill_gradient(high = "#001a33", low = "#cce6ff") 
  
  
 return(p)
}

# Test and plot three correlations ----------------------------------------

pbn_spn <- 
  pheno_cali %>%
  plot_pheno(x = pbn, y = spn)

sbn_spn <- 
  pheno_cali %>%
  plot_pheno(x = sbn, y = spn)


pbn_sbn <- 
  pheno_cali %>% 
  plot_pheno(x = pbn, y = sbn) +
  labs(caption = str_wrap("Correlation between the main traits that
                          determine panicle phenotype. A. Primary Branch 
                          Number and Spikelet Number correlates only slightly
                          mainly in wild species. B. Secondary Branch Number 
                          and Spikelet number highly correlates in cultivated
                          species, less in wild species. C. Primary and
                          Secondary Branches poorly correlates, suggesting that 
                          they are controlled by different genetic mechanisms",
                          width = 70)) +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1))


# Arrange and save plots --------------------------------------------------------------

pdf(file = "../fig/suppl-fig-correlation-pbn-spn.pdf", 
    width = 8,
    height = 9)
grid.arrange(pbn_spn,
             sbn_spn,
             pbn_sbn,
             layout_matrix = matrix(c(1,1,2,2,3,3,3), ncol = 1))
dev.off()



