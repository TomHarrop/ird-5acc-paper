library(tidyverse)
library(glue)

pc <- readRDS("data/tmp_rlog_pca_output/pc.Rds")
pcx <- readRDS("data/tmp_rlog_pca_output/pcx.Rds")
color_palette <- c("blue", "goldenrod")

# Only the first 5 PC are interesting? ------------------------------------

pc_var <- (pc$sdev/sum(pc$sdev))*100

pc_var <- tibble(var = round(pc_var,digits = 1),
                 PC_id = paste0("PC", 1:30))


pcx <- pcx %>%
  select_at(vars(accession:PC5)) %>%
  gather(PC1:PC5, key = "PC_id", value = "pc_x") %>%
  left_join(pc_var) %>%
  mutate(stage = as.character(stage),
         accession = as.character(accession)) %>%
  mutate(stage = case_when(stage == "PBM" ~ "BM",
                           stage == "SM" ~ "SM",
                           TRUE ~ stage),
         accession = case_when(accession == "japonica" ~ "O. japonica",
                               accession == "barthii" ~ "O. barthii",
                               accession == "glaberrima" ~ "O. glaberrima",
                               accession == "rufipogon" ~ "O. rufipogon",
                               accession == "indica" ~ "O. indica",
                               TRUE ~ accession)) %>%
  mutate(ID = paste0(ID, " (",
                    stage, ")"),
         accession = factor(accession,
                            levels = c("O. japonica",
                                       "O. barthii",
                                       "O. glaberrima",
                                       "O. rufipogon",
                                       "O. indica")),
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
  scale_fill_manual(values = color_palette) +
  # scale_fill_viridis_d(begin = .2, end = .8) +
  theme_bw() +
  labs(x = "Sample ID",
       y = "Loading Vector",
       fill = "Stage",
       caption = glue("Fig 2: Grafical representation of PC loadings
                       for all RNA-seq samples. The first 4 components
                       collect {sum(pc_var$var[1:4])}% of the variance
                       which is explained by differences among species.
                       The 5th component, which explains
                      {pc_var$var[5]}% of the variance, splits developmental
                      stages.") %>%
         str_wrap(width = 70)) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = .5,
                                   angle = 270),
        strip.text.x = element_text(face = "italic")) 
  

pdf("../fig/fig-02-PCA-rlog-VISUAL.pdf")
print(p)
dev.off()

