library(tidyverse)
library(readxl)

#' get_expression extracts a set of locus ids from a deseqdataset
#' object and returns them in a ggplot2 friendly format
#' 
#' Requires Tidyverse
#' 
#' Test
#' load("data/02-exp-stat.Rdata")
#' load("data/07-top-pc5-300.Rdata")
#' mapman <- read_excel("data/MAPMAN BIN-Osa_MSU_v7.xlsx")
#' 
#' dat <- get_expression(locus_ids = top_pc5,
#'                dds = dds,
#'                mapman = mapman)


get_expression <- function(locus_ids, dds, mapman = NULL) {
  
  # filter locus_ids
  excluded <- locus_ids[! locus_ids %in% rownames(dds)]
  if(length(excluded) > 0) { 
    warning(paste("These locus ids are not contained in the dds dataset:",
                  paste(excluded, collapse = ", ")))
  }
  
  locus_ids <- locus_ids[locus_ids %in% rownames(dds)]
  # first subset the DeSeqdataset, it makes everything
  # Downstream lighter
  dds <- dds[locus_ids, ]
  # print(dds)
  
  # We decided to variance stabilizing normalized counts
  dat <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  
  # Use locus id as rownames
  dat$locus_id <- rownames(dat)
  
  
  # get sample annotations
  cdata <- as.data.frame(SummarizedExperiment::colData(dds))
  cdata$sample <- rownames(cdata)
  
  # wrangle and merge
  dat <- dat %>%
    gather(colnames(dds), key = "sample", value = "normalized expression") %>%
    left_join(cdata, by = "sample") %>%
    dplyr::rename(species = accession)
  # str(dat)
  
  # merge mapman ids
  if(!is.null(mapman)) {
    mapman <- mapman %>% 
      dplyr::rename(locus_id = IDENTIFIER) %>%
      select(locus_id, DESCRIPTION) %>%
      distinct()
    # str(mapman)
    dat <- dat %>%
      left_join(mapman, by = "locus_id")
  }
  
  # scale by locus
  dat <- dat %>%
    group_by(locus_id) %>%
    mutate(by_locus = scale(`normalized expression`)) %>%
    ungroup()
  
  # scale by locus and species
  dat <- dat %>%
    group_by(species, locus_id) %>%
    mutate(by_locus_species = scale(`normalized expression`)) %>%
    ungroup()
  
  return(dat)
}


#' Load Mapman
#' 
#' requires readxl and tidyverse package
#' 
#' The Load Mapman function load and wrangles mapman files so that they can
#' be used for testing category enrichement
#' 
#' @param file a mapman file in the excel format
#' 
#' Test
#' mapman <- load_mapman(file = "data/MAPMAN BIN-Osa_MSU_v7.xlsx")

load_mapman <- function(path)
{
  mapman <- read_excel(path)
  mapman <- mapman %>% filter(! is.na(IDENTIFIER))
  
  # How many levels in the bins?
  max_l <- max(sapply(str_split(mapman$NAME, "\\."), length))
  
  # So, make it a data frame
  mapman_names <- str_split_fixed(mapman$NAME, "\\.", max_l) 
  mapman_names[nchar(mapman_names) == 0] <- NA 
  rownames(mapman_names) <- mapman$IDENTIFIER
  
  # you need some sort of cumulative string split
  # This is a very inelegant solution
  mapman_bins <- seq_along(rownames(mapman_names)) %>% #t(mapman_names)
    map(~map_chr(1:7, function(n) {
      gsub(".NA", "",
           paste(mapman_names[., 1:n], collapse = "."),
           fixed = TRUE)
    })) %>%
    do.call(rbind, .)
  mapman_bins[is.na(mapman_names)] <- NA
  mapman_bins <- data.frame(mapman_bins, stringsAsFactors = FALSE)
  colnames(mapman_bins) <- paste0("level", 1:ncol(mapman_bins))
  mapman_bins$IDENTIFIER <- mapman$IDENTIFIER
  
  # merge back
  # I would use full_join(), but IDENTIFIER are dublicated
  mapman <- bind_cols(mapman, mapman_bins)
  
  return(mapman)
}

#' `test_mapman`` tests mapman enrichement on a set of ids
#' 
#' requires tidyverse
#' 
#' @param ids a vector of locus ids
#' @param mapman a mapman dataframe, the output of `load_mapman()``
#' @param level numeric, 1 to n, the mapman level that you want to test,
#' defaults to 1
#' 
#' Test
#' load("data/07-top-pc5-300.Rdata")
#' load("data/06-mapman_bins.Rdata")
#' tst <- test_mapman(ids = top_pc5,
#'             mapman = mapman,
#'             level = 2)

test_mapman <- function(ids, mapman, level = 1)
{
  level <- paste0("level", level)
  
  # remove duplicated records in selected level
  # Otherwise statistical test overestimates entichement
  mapman_sel <- mapman %>%
    select(IDENTIFIER, level) %>%
    distinct()
  mapman_sel <- mapman_sel[complete.cases(mapman_sel), ]
  
  # Which IDENTIFIER is in ids?
  cont <- mapman_sel %>%
    select(IDENTIFIER, level) %>%
    distinct() %>%
    mutate(in_ids = ifelse(IDENTIFIER %in% ids, "yes", "no"))
  
  # Save locus_id that belong to a mapman level
  level_locus <- cont %>%
    filter(in_ids == "yes") %>%
    group_by_(level) %>%
    summarise(locus_in = paste(IDENTIFIER, collapse = " ")) %>%
    rename_(bins = level)
  
  # make it a contingency table
  cont <- table(cont[, level, drop = TRUE], cont$in_ids) %>%
    as_tibble() %>% 
    spread(key = Var2, value = n) %>%
    dplyr::rename(bins = Var1) %>%
    mutate(all_list = nrow(mapman_sel),
           in_ids = sum(cont$in_ids == "yes"))
  
  
  # Define phyper wrapper that contains "..."
  # So that it can be used in pmap with extra variables
  phyper2 <- function(q, m, n, k, ...) phyper(q, m, n, k, lower.tail = FALSE)
  
  # Test enrichment
  # inspired from
  # https://github.com/GuangchuangYu/DOSE/blob/master/R/enricher_internal.R
  cont <- cont %>%
    mutate(q = yes,  # white balls drawn
           m = no + yes, # White balls 
           n = all_list - no - yes, # Black balls
           k = in_ids) %>% # balls drawn
    # select(q, m, n, k) %>% 
    mutate(pval = pmap(., .f = phyper2, lower.tail = FALSE)) %>%
    mutate(pval = as.numeric(pval)) %>%
    # remove categories that are not detected
    filter(yes > 0) %>%
    # arrange pvalues
    arrange(pval) %>%
    left_join(level_locus, by = "bins")
  
  return(cont)
}

#' Plots normalized expression 
#' 
#' takes the output of `get_expression` 
#' and plots it with standard behaviour
#' 
#' Test
#' load("data/02-exp-stat.Rdata")
#' load("data/07-top-pc5-300.Rdata")
#' mapman <- read_excel("data/MAPMAN BIN-Osa_MSU_v7.xlsx")
#' set.seed(1)
#' dat <- get_expression(locus_ids = top_pc5,
#'                      dds = dds,
#'                      mapman = mapman) %>%
#'   filter(locus_id %in% sample(locus_id, 5))
#' p <- plot_norm_expr(dat)
#' p


plot_norm_expr <- function(dat) {
  # Set an enstablished colour palette
  color_palette <- c("blue", "goldenrod")
  
  # and a predefined species order
  dat <- dat %>% 
    mutate(species = factor(species,
                            levels = c("japonica",
                                       "barthii",
                                       "glaberrima",
                                       "rufipogon",
                                       "indica")))
  
  p <- ggplot(dat, aes(x = species,
                       y = `normalized expression`,
                       colour = stage,
                       pch = domestication)) +
    geom_point(size = 3, 
               position = position_dodge(width = .5),
               alpha = .8) +
    scale_color_manual(values = color_palette) +
    # facet_wrap(facets = c("locus_id", "DESCRIPTION"),
    #            scales = "free_y",
    #            ncol = 5,
    #            labeller = label_wrap_gen(width = 50)) +
    facet_wrap(facets = "locus_id") +
    expand_limits(y=0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotate("rect",
             xmin = 1.5, xmax = 3.5,
             ymin = -Inf, ymax = Inf,
             alpha = .1)
  return(p)
}

#' plot a keyword from gsea enrichement
#' 
#' This is a wrapper for plot_norm_expr that 
#' extracts the locus_id from a gsea enrichement element and plots them
#' 
#' test

#' test_fgsea Gets Category Enrichement on ranked gene list and Plots them
#' 
#' It is a wrapper for `fgsea` and the  from the `plotGseaTable`
#' from the `fgsea` package
#' 
#' It takes a ranking and a category as input

test_gsea <- function(rnk, mapman_list, plot_top = FALSE)
{
  gsea_res <- fgsea(pathways = mapman_list, 
                    stats = rnk,
                    minSize=15,
                    maxSize=500,
                    nperm=10000) %>% 
    arrange(padj)
  print(gsea_res$pathway[1:40])
  if(plot_top) {
    plotGseaTable(pathways = mapman_list[gsea_res$pathway[1:40]],
                  stats = sort(rnk, decreasing = T), 
                  fgseaRes = gsea_res,
                  gseaParam = 0.5)
  }
  return(gsea_res)
}
