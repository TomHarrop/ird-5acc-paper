library(tidyverse)
library(readxl)


# Wrangle Montpellier Phenotype data --------------------------------------

pheno_mnp <- read_excel("../data-raw/Phenotype_PanicleSequenced_corrected.xlsx")

# Assign descriptive names to accessions

acc <- c(B88 = "barthii",
         IR64 = "indica",
         Niponbarre = "japonica",
         Tog5681 = "glaberrima",
         W1654 = "rufipogon")

pheno_mnp <- pheno_mnp %>%
  dplyr::rename(species = `Accession Name`) %>%
  mutate(species_simple = factor(unname(acc[species]),
                                 levels = c("japonica",
                                            "barthii",
                                            "glaberrima",
                                            "rufipogon",
                                            "indica")))


# And give easier-to-type names to measurements

pheno_mnp <-  pheno_mnp %>% 
  dplyr::rename(pbn = `Pb_nb (PbN)`,
                sbn = `Sb_nb (SbN)`,
                pbl = `Pb_Length_average (PbL in cm)`,
                spn = `Sp_nb (SpN)`)

# Wrangle Cali Phenotype data ---------------------------------------------

pheno_cali <- read.table(file = "../data-raw/OsOgObOrPTRAPdata_PaperTom.txt",
                        header=T, sep="\t",
                        na.strings="NA",
                        dec=",",
                        strip.white=T,
                        stringsAsFactors = FALSE) %>% as_tibble()

# Assign simple names to  groups

acc_cali <- c(Ob = "barthii",
              Os = "sativa",
              Og = "glaberrima",
              Or = "rufipogon")

pheno_cali <- pheno_cali %>%
  mutate(Species = factor(unname(acc_cali[Origin]),
                          levels = c("barthii",
                                     "glaberrima",
                                     "rufipogon",
                                     "sativa")))

# And easier-to-type names to vairables
pheno_cali <- pheno_cali %>%
  dplyr::rename(pbn = PbN,
                sbn = SbN,
                pbl = PbL,
                spn = SpN) 


# Save Tidy Phenotype -----------------------------------------------------

save(pheno_cali, pheno_mnp, file = "../data/phenotypes.Rdata")

