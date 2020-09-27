library(dplyr)
library(erba)
library(here)
library(readxl)

# Supplementary Figure 1

##########################################################################
# A) KOs Transcription Factors #####

## load data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                sheet = 1,
                                col_names = TRUE,
                                skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "general",
                  column_total = total_kos_tf$total,
                  column_orfs = total_kos_tf$`ORFs(X100)`,
                  filename = here::here("figures/S1_kos_tf_lm.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome")

### obtain lm Coefficients in general
summary(lm(total_kos_tf$total ~ total_kos_tf$`ORFs(X100)`) )
cor(total_kos_tf$total , total_kos_tf$`ORFs(X100)`)

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
                  column_total = total_kos_sigma$total,
                  column_orfs = total_kos_sigma$`ORFs(X100)`,
                  filename =  here::here("figures/S1_kos_sf_lm.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  ymax = 120)

## obtain lm Coefficients in general
summary(lm(total_kos_sigma$total ~ total_kos_sigma$`ORFs(X100)`) )
cor(total_kos_sigma$total , total_kos_sigma$`ORFs(X100)`)
###############################################################################
# C) KOs Transcription Factors per phylum #####

## load data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                   sheet = 1,
                                   col_names = TRUE,
                                   skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_kos_tf, type = "groups",
                  column_total = total_kos_tf$total,
                  column_orfs = total_kos_tf$`ORFs(X100)`,
                  filename = here::here("figures/S1_kos_tf_lm_color.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome")

### obtain lm Coefficients per phylum
erba::get_slopePerPhylum(total_kos_tf, x = "ORFs(X100)", y = "total")
erba::get_correlation(total_kos_tf, x = "ORFs(X100)", y = "total")


get_R2_perPhylum <- function(x) {
    y <- dplyr::filter(total_kos_tf, phylum  == x)
    summary(lm(y$total ~ y$`ORFs(X100)`))
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
                  column_total = total_kos_sigma$total,
                  column_orfs = total_kos_sigma$`ORFs(X100)`,
                  filename =  "figures/S1_kos_sf_lm_color.tiff",
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  ymax = 120)

## obtain lm Coefficients per phylum
erba::get_slopePerPhylum(total_kos_sigma, x = "ORFs(X100)", y = "total")
erba::get_correlation(total_kos_sigma, x = "ORFs(X100)", y = "total")

get_R2_perPhylum <- function(x) {
    y <- dplyr::filter(total_kos_sigma, phylum  == x)
    summary(lm(y$total ~ y$`ORFs(X100)`))
}

sapply(sort(unique(total_kos_sigma$phylum)), get_R2_perPhylum)

##########################################################################
