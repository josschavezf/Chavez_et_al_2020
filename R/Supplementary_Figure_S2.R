library(dplyr)
library(erba)
library(here)
library(readxl)

# Supplementary Figure 2

##########################################################################
# A) KOs Transcription Factors #####

## load data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                sheet = 1,
                                col_names = TRUE,
                                skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "general",
                  x = total_kos_tf$CDS/100,
                  y = total_kos_tf$total,
                  filename = here::here("figures/S2_kos_tf_lm.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome",
                  xlab = "CDS (x 100)")

### obtain lm Coefficients in general
total_kos_tf$CDS <- total_kos_tf$CDS/100
summary(lm(total_kos_tf$total ~ total_kos_tf$CDS) )
cor(total_kos_tf$total , total_kos_tf$CDS) %>% round(2)

####################################################################
# B) KOs Sigma Factors ####

## load data
total_kos_sigma <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                   sheet = 2,
                                   col_names = TRUE,
                                   skip = 2)
archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_kos_sigma <- total_kos_sigma %>%
    filter(!phylum %in% archaeas)

## plot distribution of sigma factors versus genome size
erba::plot_points(total_kos_sigma, type =  "general",
                  x = total_kos_sigma$CDS/100,
                  y = total_kos_sigma$total,
                  filename =  here::here("figures/S2_kos_sf_lm.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  xlab = "CDS (x 100)",ymax = 120)

## obtain lm Coefficients in general
total_kos_sigma$CDS <- total_kos_sigma$CDS/100
summary(lm(total_kos_sigma$total ~ total_kos_sigma$CDS) )
cor(total_kos_sigma$total , total_kos_sigma$CDS) %>% round(2)

###############################################################################
# C) KOs Transcription Factors per phylum #####

## load data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                   sheet = 1,
                                   col_names = TRUE,
                                   skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "groups",
                  x = total_kos_tf$CDS/100,
                  y = total_kos_tf$total,
                  filename = here::here("figures/S2_kos_tf_lm_color.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome",
                  xlab = "CDS (x 100)")

### obtain lm Coefficients per phylum
total_kos_tf$CDS <- total_kos_tf$CDS/100
erba::get_slopePerPhylum(total_kos_tf, x = "CDS", y = "total")
erba::get_correlation(total_kos_tf, x = "CDS", y = "total")


get_R2_perPhylum <- function(x) {
    y <- dplyr::filter(total_kos_tf, phylum  == x)
    summary(lm(y$total ~ y$CDS))
}

sapply(sort(unique(total_kos_tf$phylum)), get_R2_perPhylum)

####################################################################
# D) KOs Sigma Factors per phylum ####

## load data
total_kos_sigma <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)
archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_kos_sigma <- total_kos_sigma %>%
  filter(!phylum %in% archaeas)

## plot distribution of sigma factors versus genome size
erba::plot_points(total_kos_sigma, type =  "groups",
                  x = total_kos_sigma$CDS/100,
                  y = total_kos_sigma$total,
                  filename =  "figures/S2_kos_sf_lm_color.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  xlab = "CDS (x 100)", ymax = 120)

## obtain lm Coefficients per phylum
total_kos_sigma$CDS <- total_kos_sigma$CDS/100
erba::get_slopePerPhylum(total_kos_sigma, x = "CDS", y = "total")
erba::get_correlation(total_kos_sigma, x = "CDS", y = "total")

get_R2_perPhylum <- function(x) {
    y <- dplyr::filter(total_kos_sigma, phylum  == x)
    summary(lm(y$total ~ y$CDS))
}

sapply(sort(unique(total_kos_sigma$phylum)), get_R2_perPhylum)

##########################################################################
