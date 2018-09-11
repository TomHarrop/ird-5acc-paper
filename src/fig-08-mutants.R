library(tidyverse)
library(readxl)
library(ggpubr)
library(gridExtra)

dat <- read_excel(path = "../data-raw/Mutant_AP2_PhenotypingData.xlsx")
gene_names <- c(LOC_Os07g03250 = "plt8",
                WT = "WT",
                LOC_Os05g32270 = "erf142",
                LOC_Os06g03710 = "smo2 ???",
                LOC_Os08g31580 = "erf48")


dat <- dat %>% 
  mutate_at(vars(RL:SpN), as.numeric) %>% 
  mutate(mutant_gene = gene_names[`Target Gene`]) %>%
  gather(RL:SpN, key = "measure", value = "value")


# Plot every mutant -------------------------------------------------------

# pdf("mutants/mutants-boxplot.pdf", height = 14, width = 5)
ggplot(dat, aes(x = ID, y = value)) +
  geom_boxplot(varwidth = T) +
  facet_grid(measure ~ Accession, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()


# Test --------------------------------------------------------------------

tst <- dat %>%
  filter(Accession != "Kitake") %>%
  arrange(desc(mutant_gene)) %>%
  mutate(mutant_gene = as_factor(mutant_gene)) %>%
  split(paste(.$measure, .$Accession)) %>%
  map(~lm(value ~ mutant_gene, data = .)) %>%
  map(summary)

# save(tst, file = "mutants/tests.Rdata")


# Test with r4ds functions ------------------------------------------------

test_pheno <- function(df) {
  lm(value ~ mutant_gene, data = df)
}

# Not sure how to extract coefficients systematically

tst <- dat %>%
  filter(Accession != "Kitake") %>%
  arrange(desc(mutant_gene)) %>%
  mutate(mutant_gene = as_factor(mutant_gene)) %>%
  group_by(measure, Accession) %>%
  nest() %>%
  mutate(model = map(data, test_pheno)) %>%
  mutate(model_summ = map(model, summary))
# mutate(pvals = map(model_summ, .f = function(i) i$a$coefficients[, "Pr(>|t|)"]))
# mutate(pval)


# How to add pvalue to a plot? (ggpubr)--------------------------------------

tst <- dat %>%
  split(.$Accession)

tst <- dat %>%
  filter(Accession != "Kitake") %>%
  arrange(desc(mutant_gene)) %>%
  mutate(mutant_gene = as_factor(mutant_gene)) %>%
  split(paste(.$Accession)) %>%
  map(~split(., .$measure))

p <- list()
p$Nipponbare <- map(tst$Niponbarre, ~ggboxplot(., x = "mutant_gene", y = "value",
                                               ylab = unique(.$measure),
                                               title = unique(.$measure),
                                               fill = "cornsilk") +
                      stat_compare_means(method = "t.test",
                                         comparisons = list(c(1,2), c(2,3), c(1,3))))

p$Kinmaze <- map(tst$Kinmaze, ~ggboxplot(., x = "mutant_gene", y = "value",
                                         ylab = unique(.$measure),
                                         title = unique(.$measure),
                                         fill = "cornsilk") +
                   stat_compare_means(method = "t.test",
                                      comparisons = list(c(1,2))))

p$Illmibyeo <- map(tst$Illmibyeo, ~ggboxplot(., x = "mutant_gene", y = "value",
                                             ylab = unique(.$measure),
                                             title = unique(.$measure),
                                             fill = "cornsilk") +
                     stat_compare_means(method = "t.test",
                                        comparisons = list(c(1,2))))

arrange_plots <- function(name) {
  header <- text_grob(name,
                      face = "bold",
                      size = 32)
  p1 <- ggarrange(plotlist = p[[name]],
                  ncol = length(p[[name]])/3,
                  nrow = length(p[[name]])/3)
  ggarrange(header, p1, nrow = 2, 
            heights = c(1, 9))
}

# pdf(file = "mutants/mutants_with_stats.pdf",
    # width = 12,
    # height = 14)
map(names(p), arrange_plots)
# dev.off()

# Check prediction --------------------------------------------------------

load("data/phenotype-fit.Rdata")

dat <- read_excel(path = "data/Mutant_AP2_PhenotypingData.xlsx")

dat <- dat %>% 
  mutate_at(vars(RL:SpN), as.numeric) %>%
  rename(pbn = PbN,
         sbn = SbN,
         spn = SpN) %>%
  mutate(spn_pred = predict(fit, .))

ggplot(dat, aes(x = spn_pred,
                y = spn)) +
  geom_point(alpha = .3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_grid(`Target Gene` ~ Accession) +
  theme_bw()
# dat <- dat %>%
#   split(.$measure)
# pdf("figures/mutants-boxplot.pdf", height = 4)
# names(dat) %>% 
#   map(~ggplot(dat[[.]], aes(x = `Target Gene`, y = value)) +
#         geom_boxplot() +
#         # geom_point(size = 2, alpha = .2) +
#         facet_wrap(facets = "Accession",
#                    nrow = 1) +
#         ggtitle(.) +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)))
# dev.off()