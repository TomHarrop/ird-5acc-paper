library(tidyverse)
library(glue)

load("../data/rlog-pca.Rdata")
# color_palette <- c("blue", "goldenrod")

# Only the first 5 PC are interesting? ------------------------------------

pc_var <- (pc$sdev/sum(pc$sdev))*100

pc_var <- tibble(var = round(pc_var,digits = 1),
                 PC_id = paste0("PC", 1:30))


pcx <- pcx %>%
  select_at(vars(accession:PC5)) %>%
  gather(PC1:PC5, key = "PC_id", value = "pc_x") %>%
  left_join(pc_var) %>%
  mutate(ID = paste0(ID, " (",
                    stage, ")"),
         accession = factor(accession,
                            levels = c("japonica",
                                       "barthii",
                                       "glaberrima",
                                       "rufipogon",
                                       "indica")),
         PC_id = paste0(PC_id, "~", var, "%"))


# Plot and save -----------------------------------------------------------

p <- ggplot(pcx,
       aes(x = ID,
           fill = stage,
           # colour = stage,
           y = pc_x)) +
  geom_col() +
  # geom_point(size = 4) +
  geom_hline(yintercept = 0) +
  facet_grid(PC_id ~ accession, scales = "free") +
  # scale_fill_manual(values = color_palette) +
  scale_fill_viridis_d(begin = .2, end = .8) +
  theme_bw() +
  labs(x = "Sample ID",
       y = "Loading",
       fill = "Stage",
       caption = glue("Fig 2: Grafical representation of PC loadings
                       for all RNA-seq samples. The first 4 components
                       explain {sum(pc_var$var[1:4])}% of the variance
                       which might be determined by species specific 
                      differences. Only the 5th component, which explains
                      {pc_var$var[5]}% of the splits stages.") %>%
         str_wrap(width = 70)) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = .5,
                                   angle = 270)) 
  
  
pdf("../fig/fig-02-PCA-rlog-VISUAL.pdf")
print(p)
dev.off()

