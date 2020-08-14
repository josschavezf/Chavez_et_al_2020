library(erba)
library(tidyverse)
library(here)

# Supplementary Figure 1
###############################################################################
# KOs Transcription Factors #####

## load data
total_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")
total_kos_tf <-

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "general",
                  filename = "figures/kos_tf_factor_lm.tiff",
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome")

### obtain lm Coefficients in general
lm(data = data_kos_tf, formula = total~ORFs)
cor(total_kos_tf$total, total_kos_tf$ORFs) %>% round(2)

####################################################################
# KOs Sigma Factors ####

## load data
total_kos_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_9.txt")
archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_kos_sigma <- total_kos_sigma %>%
  filter(!phylum %in% archaeas)

## plot distribution of sigma factors versus genome size
erba::plot_points(total_kos_sigma, type =  "general",
                  filename =  "figures/kos_sigma_factor_lm.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  ymax = 120)

## obtain lm Coefficients in general
lm(data = total_kos_sigma, formula = total~ORFs)
cor(total_kos_sigma$total, total_kos_sigma$ORFs) %>% round(2)

##########################################################################
