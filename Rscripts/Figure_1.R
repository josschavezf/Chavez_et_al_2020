library(erba)
library(tidyverse)

##########################################################################
# COGs Transcription Factors #####

## load data
total_cogs_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_6.txt")

## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs_tf, type = "general",
                  filename = "figures/cogs_tf_factor_lm.tiff",
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome", ymax = 800)

### obtain lm Coefficients in general
lm(data = total_cogs_tf, formula = total~ORFs)
cor(total_cogs_tf$total, total_cogs_tf$ORFs) %>% round(2)
####################################################################
# COGs Sigma Factors ####

## load data
total_cogs_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_8.txt")

archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_cogs_sigma <- total_cogs_sigma %>%
  filter(!phylum %in% archaeas)

## plot distribution of sigma factors versus genome size
erba::plot_points(total_cogs_sigma, type =  "general",
                  filename =  "figures/cogs_sigma_factor_lm.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 150)

## obtain lm Coefficients in general
lm(data = total_cogs_sigma, formula = total~ORFs)
cor(total_cogs_sigma$total, total_cogs_sigma$ORFs) %>% round(2)

##########################################################################
# KOs Transcription Factors #####

## load data
total_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")

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
# Riboswitch #####

## load data
data_riboswitch <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_10.txt")

## plot distribution of riboswitches versus genome size
erba::plot_exception(data_riboswitch,
                     filename = "figures/riboswitch_lm.tiff",
                     title = "Transcriptional Riboswitches",
                     ylab = "Riboswitches per genome",
                     ymax = 80,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",], formula = total~ORFs)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",], formula = total~ORFs)

data_riboswitch  %>%
    filter(phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)

data_riboswitch  %>%
    filter(!phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)
