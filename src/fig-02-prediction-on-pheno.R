library(tidyverse)
library(ggfortify)
library(gridExtra)

load("../data/phenotypes.Rdata")

pbn_spn <- pheno_cali %>% 
  ggplot(aes(x = pbn, y = spn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 

sbn_spn <- pheno_cali %>% 
  ggplot(aes(x = sbn, y = spn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") 

pbn_sbn <- pheno_cali %>% 
  ggplot(aes(x = pbn, y = sbn)) +
  # geom_point(alpha = .2) +
  geom_hex(bins = 20) +
  facet_wrap(facets = "Species",
             nrow = 1) +
  theme_bw() +
  labs(caption = str_wrap("Correlation between the main traits that
                          determine panicle phenotype. A. Primary Branch 
                          Number and Spikelet Number correlates only slightly
                          mainly in wild species. B. Secondary Branch Number 
                          and Spikelet number highly correlates in cultivated
                          species, less in wild species. C. Primary and
                          Secondary Branches poorly correlates, suggesting that 
                          they are controlled by different genetic mechanisms",
                          width = 70))+
  scale_fill_gradient(high = "#001a33", low = "#cce6ff") +
  theme(plot.caption = element_text(size = 15, hjust = 0, lineheight = 1))

pdf(file = "../fig/fig-02-descript-pheno-stat.pdf", 
    width = 8,
    height = 9)
grid.arrange(pbn_spn,
             sbn_spn,
             pbn_sbn,
             layout_matrix = matrix(c(1,1 ,2,2,3,3,3), ncol = 1))
dev.off()
# # Old ----------------------------------------------------

# 
# pheno_cali <- pheno_cali %>%
#   mutate(bn = pbn + sbn + TbN) %>%
#   mutate(naxil = (spn - pbn)/(pbn + sbn +TbN))
# 
# # 01 - Plot axillary meristems per branch --------------------------------
# 
# plot_naxil <- pheno_cali %>%
#   ggplot(aes(x = spn,
#              y = naxil,
#              colour = Type,
#              pch = Species)) + 
#   geom_point(alpha = .5) +
#   ylab("Number of Axillary Meristems") + 
#   xlab("Spikelet Number") +
#   # facet_wrap(facets = "Species") +
#   theme_bw()
# 
# # 02 - Domesticated species produce more primary branches -----------------
# 
# t.test(pbn ~ Type, data = pheno_cali) 
# 
# # 03 - Secondary Branches correlates With Spikelet Number------------------
# 
# plot_corr <- ggplot(pheno_cali, aes(y = sbn,
#                        x  = spn,
#                        colour = Type,
#                        pch = Species)) + 
#   geom_point(alpha = .5) +
#   # geom_hex(bins = 30) +
#   # facet_wrap(facets =  "Species") +
#   theme_bw() +
#   xlab("Spikelet Number") + 
#   ylab("Secondary Branch Number")
#   # scale_fill_gradient(high = "#001a33", low = "#cce6ff") 
# 
# # 04 - Fit many models ----------------------------------------------------
# 
# fit_cali <- pheno_cali %>% lm(spn ~ pbn + sbn + TbN + 0,
#                                   data = .)
# # summary(fit_cali)
# 
# # fit_cali_int <- pheno_cali %>% lm(spn ~ pbn + sbn + pbn*sbn + 0,
# #                    data = .)
# # summary(fit_cali_int)
# # 
# # fit_cali_all <- pheno_cali %>% lm(spn ~ pbn + sbn + pbn*sbn + TbN + 0,
# #                                   data = .)
# # summary(fit_cali_all)
# # 
# # 
# # fit_cali_bn <- pheno_cali %>% lm(spn ~ bn + 0,
# #                                  data = .)
# # 
# # fit_sbn <- pheno_cali %>% lm(spn ~ sbn + 0,
# #                              data = .)
# 
# 
# # 05 - Predict and Plot --------------------------------------------------------
# 
# pheno_mnp <- pheno_mnp %>%
#   dplyr::rename(TbN = "Tb_nb (TbN)") %>%
#   mutate(bn = pbn + sbn + TbN)
# 
# pheno_mnp <- pheno_mnp %>%
#   mutate(spn_pred = predict(fit_cali,
#                             newdata = pheno_mnp)) %>%
#   mutate(Type = ifelse(species_simple %in% c("barthii", "rufipogon"),
#                        "Wild","Cultivated"))
# 
# 
# pred_plot <- pheno_mnp %>% 
#   ggplot(aes(x = spn,
#              y = spn_pred, 
#              colour = Type,
#              pch = Species
#   )) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   xlab("Spikelet Number - Observed") + 
#   ylab("Spikelet Number - Predicted from Branch Number")
#   # facet_wrap(facets = "Species", nrow = 1) +
#   theme_bw()
# 
# 
# # Caption -----------------------------------------------------------------
# 
# capt <- str_wrap("Rice lodges more Spikelet on its panicle by
#                  increasing branch number and not by increasing 
#                  spikelets per branch. A. Spikelet per Branch in high 
#                  yield accession is comparable to the low yield wild
#                  accessions. B. Spikelet Number correlates with secondary branch
#                  number, in most cultivated species. Low yielding species might 
#                  produce outliers that are better explained bt primary branch number
#                  C. A linear model with SpN as outcome and PbN + Sbn + TbN
#                  as predictors is able to predict spikelet number is an independent
#                  dataset", width = 80)  
# 
# # Put together plots and save ---------------------------------------------
# 
# pdf(file = "../fig/fig-02-prediction-on-pheno.pdf", 
#     width = 12,
#     height = 3)
# grid.arrange(plot_naxil,
#              plot_corr,
#              pred_plot,
#              layout_matrix = matrix(c(1,1,1,2,2,2,3,3,3,3), nrow = 1))
# dev.off()
# 
# # Old stuff ---------------------------------------------------------------
# 
# 
# # pheno_mnp %>% 
# #   mutate(spn_pred = predict(fit_cali_int, newdata = pheno_mnp)) %>%
# #   ggplot(aes(x = spn, y = spn_pred)) +
# #   geom_point() +
# #   geom_abline(slope = 1, intercept = 0) +
# #   facet_wrap(facets = "Species") + 
# #   theme_bw()
# # 
# # 
# # pheno_mnp %>% 
# #   mutate(spn_pred = predict(fit_cali_bn, newdata = pheno_mnp)) %>%
# #   ggplot(aes(x = spn, y = spn_pred)) +
# #   geom_point() +
# #   geom_abline(slope = 1, intercept = 0) +
# #   facet_wrap(facets = "Species")
# # 
# # 
# # pheno_mnp %>% 
# #   mutate(spn_pred = predict(fit_sbn, newdata = pheno_mnp)) %>%
# #   ggplot(aes(x = spn, y = spn_pred)) +
# #   geom_point() +
# #   geom_abline(slope = 1, intercept = 0) +
# #   facet_wrap(facets = "Species") + 
# #   theme_bw()
# # 
# # 
# # pheno_mnp%>% lm(spn ~ pbn + sbn + 0,
# #                 data = .) %>%
# #   summary()
# # 
