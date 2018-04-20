library(tidyverse)
library(readxl)
library(gridExtra)
library(leaps)

load("../data/phenotypes.Rdata")


# Search for the best model on mnp data -----------------------------------


reg_mnp <- regsubsets(spn ~ . + 0, pheno_mnp %>%
                        select_at(vars(`Rachis_length (RL)`:spn)) %>%
                        # no tertiary branches detected in this dataset
                        select(-`Tb_nb (TbN)`))

reg_mnp <- summary(reg_mnp)
# plot (reg_mnp$cp, xlab = "Number of Variables" , ylab ="Cp" ,
#       type = "l")
which.min(reg_mnp$cp)
reg_mnp$which[which.min(reg_mnp$cp), ]

# Fit model with pbn sbn and pbl ------------------------------------------

fit <- lm(spn ~ pbn + sbn  + 0,
          data = pheno_mnp)

pheno_cali <- pheno_cali %>%
  mutate(prediction_montpellier = predict(object = fit,
                                          newdata = pheno_cali))

svg("../fig/fig-spn-prediction.svg",
    height = 5.8,
    width = 5.5)
ggplot(pheno_cali, aes(x = spn,
                       y = prediction_montpellier)) + 
  geom_point(size = 2, alpha = .1) +
  # geom_hex() +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(facets = "Species") +
  xlab("Observed") +
  ylab("Predicted") +
  ggtitle("Spikelet Number is predicted by primary and\n secondary branch number") +
  theme_bw() 
dev.off()

