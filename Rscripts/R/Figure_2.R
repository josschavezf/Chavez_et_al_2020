library(dplyr)
library(erba)
library(here)
library(readxl)

# read TFs KO description
kos_description <- readxl::read_excel(here::here("data/Table_S1.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)
kos_repressor <- kos_description %>%
  filter(Description %in%
           kos_description$Description[grep("repressor", kos_description$Description)]) %>%
  select(KO)

kos_non_repressor <- kos_description %>%
  filter(!KO %in% kos_repressor$KO ) %>%
  select(KO)

## load TF counts data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                   sheet = 1,
                                   col_names = TRUE,
                                   skip = 2)

## summarise repressor and non_repressor data

tf_repressor_non_repressor <- total_kos_tf %>%
  mutate(repressor = rowSums(total_kos_tf[all_of(kos_repressor$KO)] ),
         non_repressor = rowSums(total_kos_tf[all_of(kos_non_repressor$KO)]) ) %>%
  select("organism",  "ORFs(X100)", "repressor", "non_repressor")
colnames(tf_repressor_non_repressor)[2] <- "ORFs"

## plot data
erba::plot_repressor(tf_repressor_non_repressor,
                     filename = here::here("figures/tf_repressor_non_repressor.tiff"),
                     title = "Transcription factors per genome",
                     ylab = "Transcription factors",
                     ymax = 200)

## get linear correlation coefficients
lm(data = tf_repressor_non_repressor, formula = repressor ~ ORFs)
lm(data = tf_repressor_non_repressor, formula = non_repressor ~ ORFs)

round(cor(tf_repressor_non_repressor$repressor, tf_repressor_non_repressor$ORFs), 2)
round(cor(tf_repressor_non_repressor$non_repressor, tf_repressor_non_repressor$ORFs), 2)
