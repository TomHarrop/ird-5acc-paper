library(tidyverse)
library(ggfortify)
library(gridExtra)

load("../data/phenotypes.Rdata")

pheno_cali <- pheno_cali %>%
  mutate(bn = pbn + sbn + TbN)



# 01 - Domesticated species produce more primary branches -----------------

t.test(pbn ~ Type, data = pheno_cali) 

# 02 - A model with 

fit_cali <- pheno_cali %>% lm(spn ~ pbn + sbn + TbN + 0,
                                  data = .)
summary(fit_cali)

fit_cali_int <- pheno_cali %>% lm(spn ~ pbn + sbn + pbn*sbn + 0,
                   data = .)
summary(fit_cali_int)

fit_cali_all <- pheno_cali %>% lm(spn ~ pbn + sbn + pbn*sbn + TbN + 0,
                                  data = .)
summary(fit_cali_all)


fit_cali_bn <- pheno_cali %>% lm(spn ~ bn + 0,
                                 data = .)

fit_sbn <- pheno_cali %>% lm(spn ~ sbn + 0,
                             data = .)



pheno_mnp <- pheno_mnp %>%
  dplyr::rename(TbN = "Tb_nb (TbN)") %>%
  mutate(bn = pbn + sbn + TbN)

pheno_mnp %>% 
  mutate(spn_pred = predict(fit_cali, newdata = pheno_mnp)) %>%
  ggplot(aes(x = spn, y = spn_pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(facets = "Species") + 
  theme_bw()
pheno_mnp %>% 
  mutate(spn_pred = predict(fit_cali_int, newdata = pheno_mnp)) %>%
  ggplot(aes(x = spn, y = spn_pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(facets = "Species") + 
  theme_bw()


pheno_mnp %>% 
  mutate(spn_pred = predict(fit_cali_bn, newdata = pheno_mnp)) %>%
  ggplot(aes(x = spn, y = spn_pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(facets = "Species")


pheno_mnp %>% 
  mutate(spn_pred = predict(fit_sbn, newdata = pheno_mnp)) %>%
  ggplot(aes(x = spn, y = spn_pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(facets = "Species") + 
  theme_bw()


pheno_mnp%>% lm(spn ~ pbn + sbn + 0,
                data = .) %>%
  summary()

