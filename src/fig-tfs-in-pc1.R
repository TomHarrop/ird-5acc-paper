library(tidyverse)
# library(DESeq2)
source("helper-functions.R")
library(fgsea)

load("../data/all-sep-deseq.Rdata")
dds <- readRDS("../data-raw/dds.Rds")
tf_fam <- readRDS("../data-raw/tfdb_os.Rds") %>%
  dplyr::rename(locus_id = "Protein.ID") %>%
  filter(locus_id %in% rownames(dds)) %>%
  filter(!duplicated(locus_id))

pcx <- pc_spc$x %>% 
  as.data.frame(.) %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  left_join(tf_fam) %>%
  arrange(desc(PC3)) %>%
  mutate(rank_pc3 = 1:nrow(.)) %>%
  arrange(desc(PC2)) %>%
  mutate(rank_pc2 = 1:nrow(.)) %>%
  arrange(desc(PC1)) %>%
  mutate(rank_pc1 = 1:nrow(.)) %>%
  mutate(Family = ifelse(is.na(Family), "none", Family))

pcx <- pcx %>% left_join(test_gsea(set_names(pcx$PC1, 
                                       nm = pcx$locus_id),
                             mapman = split(pcx$locus_id, pcx$Family)) %>%
                     select(pathway, padj) %>%
                     dplyr::rename(Family = "pathway",
                                   pc1_padj = "padj")) %>%
  left_join(test_gsea(set_names(pcx$PC2, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
              select(pathway, padj) %>%
              dplyr::rename(Family = "pathway",
                            pc2_padj = "padj")) %>%
  left_join(test_gsea(set_names(pcx$PC3, 
                                nm = pcx$locus_id),
                      mapman = split(pcx$locus_id, pcx$Family)) %>%
              select(pathway, padj) %>%
              dplyr::rename(Family = "pathway",
                            pc3_padj = "padj"))

svg(filename = "../fig/fig-tfs-in-pc1.svg",
    height = 2)
ggplot(pcx %>% 
         filter(Category == "TF"),
         # filter(Family == "MADS"),
       aes(x = rank_pc1,
           y = PC1)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC1), lwd = .05) + 
  geom_hline(yintercept = 0, lwd = .05) +
  geom_rug(aes(x = rank_pc1, y = NULL), alpha = .1) +
  theme_bw()
dev.off()  

pcx_tf <- pcx%>%
  filter(Category == "TF") %>%
  group_by(Family) %>%
  add_tally() %>%
  ungroup() %>%
  filter(n > 5)

pdf("../fig/fig-tfs-in-pc1.pdf",
    height = 50, width = 6)
ggplot(pcx_tf,
         # filter(Family %in% c("AP2-EREBP", "MADS", "WRKY")),
       aes(x = rank_pc1,
           y = PC1)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC1), lwd = .5) + 
  geom_hline(yintercept = 0, lwd = .05) +
  geom_rug(aes(x = rank_pc1, y = NULL), alpha = .5) +
  facet_wrap(facets = c( "pc1_padj", "Family"), ncol = 1) +
  theme_bw()
dev.off()

pdf("../fig/fig-tfs-in-pc2.pdf",
    height = 50, width = 6)
ggplot(pcx_tf,
       # filter(Family %in% c("AP2-EREBP", "MADS", "WRKY")),
       aes(x = rank_pc2,
           y = PC2)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC2), lwd = .5) + 
  geom_hline(yintercept = 0, lwd = .05) +
  geom_rug(aes(x = rank_pc2, y = NULL), alpha = .5) +
  facet_wrap(facets = c( "pc2_padj", "Family"), ncol = 1) +
  theme_bw()
dev.off()

pdf("../fig/fig-tfs-in-pc3.pdf",
    height = 50, width = 6)
ggplot(pcx_tf,
       # filter(Family %in% c("AP2-EREBP", "MADS", "WRKY")),
       aes(x = rank_pc3,
           y = PC3)) + 
  # geom_point(alpha = .01) + 
  geom_linerange(aes(ymin = 0, ymax = PC3), lwd = .5) + 
  geom_hline(yintercept = 0, lwd = .05) +
  geom_rug(aes(x = rank_pc3, y = NULL), alpha = .5) +
  facet_wrap(facets = c( "pc3_padj", "Family"), ncol = 1) +
  theme_bw()
dev.off()

