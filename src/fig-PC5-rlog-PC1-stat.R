library(tidyverse)
load("../data/all-sep-deseq.Rdata")
load("../data/rlog-pca.Rdata")


pc_spc <- pc_spc$x %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(locus_id = "rowname") %>%
  select(locus_id, PC1)

to_plot <- pcro %>%
  select(PC5, locus_id) %>%
  inner_join(pc_spc)


pdf(file = "../fig/fig-PC5-rlog-PC1-stat.pdf")
ggplot(to_plot,
       aes(x = PC5,
           y = PC1)) + 
  geom_hex(bins = 80) +
  theme_bw() +
  xlab("PC5 - rlog") +
  ylab("PC1 - stat") +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff")
dev.off()
