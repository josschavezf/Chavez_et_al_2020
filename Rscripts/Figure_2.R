library(erba)
library(dplyr)
##########################################################################

# COGs Transcription Factors #####

## load data
total_cogs_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_6.txt")

## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs_tf, type = "groups",
                  filename = "figures/cogs_tf_factor_lm_color.tiff",
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome", ymax = 800)

## obtain lm Coefficients per group
get_correlation(total_cogs_tf)
get_linear_coefficients(total_cogs_tf)

####################################################################
# COGs Sigma Factors ####

## load data
total_cogs_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_8.txt")
archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_cogs_sigma <- total_cogs_sigma %>%
  filter(!phylum %in% archaeas)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs_sigma, type =  "groups",
                  filename  =  "figures/cogs_sigma_factor_lm_color.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 150)

### obtain lm Coefficients per group
erba::get_correlation(total_cogs_sigma)
erba::get_linear_coefficients(total_cogs_sigma)

##########################################################################

# KOs Transcription Factors #####

## load data
total_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "groups",
                  filename = "figures/kos_tf_factor_lm_color.tiff",
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome")

## obtain lm Coefficients per group
get_correlation(total_kos_tf)
get_linear_coefficients(total_kos_tf)

####################################################################
# KOs Sigma Factors ####

## load data
total_kos_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_9.txt")
archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_kos_sigma <- total_kos_sigma %>%
  filter(!phylum %in% archaeas)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_sigma, type =  "groups",
                  filename  =  "figures/kos_sigma_factor_lm_color.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 120)

## obtain lm Coefficients per group
erba::get_correlation(total_kos_sigma)
erba::get_linear_coefficients(total_kos_sigma)

##########################################################################
# Riboswitch #####

## load data
data_riboswitch <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_10.txt")

## plot distribution of riboswitches versus genome size
zero_groups <- c("Chlamydiae", "Crenarchaeota")
data_riboswitch <- data_riboswitch %>%
  filter(!phylum %in% zero_groups)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  filename = "figures/riboswitch_lm_color.tiff",
                  title = "Transcriptional Riboswitches",
                  ylab = "Riboswitches per genome",
                  ymax = 80)

### obtain lm Coefficients per group
erba::get_correlation(data_riboswitch)
erba::get_linear_coefficients(data_riboswitch)
