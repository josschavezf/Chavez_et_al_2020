library(dplyr)
library(erba)
library(here)
library(readxl)

##########################################################################
# A) COGs Transcription Factors #####

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                            col_names = TRUE,
                            skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs, type = "general",
                  x = total_cogs$CDS/100,
                  y = total_cogs$`Transcription factors`,
                  filename = here::here("figures/F1_cogs_tf_lm.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome",
                  xlab = "CDS (x 100)",
                  ymax = 800)

### obtain lm Coefficients in general
total_cogs$CDS <- total_cogs$CDS/100
lm(total_cogs$`Transcription factors` ~ total_cogs$CDS)
cor(total_cogs$`Transcription factors`, total_cogs$CDS) %>% round(2)

summary(lm(total_cogs$`Transcription factors` ~ total_cogs$CDS))

####################################################################
# B) COGs Sigma Factors ####

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                 col_names = TRUE,
                                 skip = 2)

archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_cogs <- total_cogs %>%
  filter(!phylum %in% archaeas)

## plot distribution of sigma factors versus genome size
erba::plot_points(total_cogs, type =  "general",
                  x = total_cogs$CDS/100,
                  y = total_cogs$`Sigma factors`,
                  filename =  here::here("figures/F1_cogs_sf_lm.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  xlab = "CDS (x 100)",ymax = 150)

## obtain lm Coefficients in general
lm(total_cogs$`Sigma factors` ~ total_cogs$CDS )
cor(total_cogs$`Sigma factors` , total_cogs$CDS )

summary(lm(total_cogs$`Sigma factors` ~ total_cogs$CDS))
##########################################################################
# C) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

zero_riboswitches <- c("Chlamydiae", "Crenarchaeota")
data_riboswitch <- data_riboswitch %>% dplyr::filter(!phylum %in% zero_riboswitches)
data_riboswitch$CDS <- data_riboswitch$CDS/100

erba::plot_exception(data_riboswitch,
                     filename = here::here("figures/F1_riboswitch_lm.tiff"),
                     title = "Transcriptional Riboswitches",
                     ylab = "Riboswitches per genome",
                     ymax = 20,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~CDS)[[1]][2] %>% round(2)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~CDS)[[1]][2] %>% round(2)

filter(data_riboswitch, phylum == "Firmicutes") %>%
  summarise(cor = cor(total, CDS) )
filter(data_riboswitch, !phylum == "Firmicutes") %>%
  summarise(cor = cor(total, CDS) )

summary(lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
           formula = total~ORFs))
summary(lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
           formula = total~ORFs))
##########################################################################
# D) COGs Transcription Factors #####

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                 col_names = TRUE,
                                 skip = 2)

## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs, type = "groups",
                  x = total_cogs$CDS/100,
                  y = total_cogs$`Transcription factors`,
                  filename = here::here("figures/F1_cogs_tf_lm_color.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome",
                  xlab = "CDS (x 100)", ymax = 800)

## obtain lm Coefficients per group
total_cogs$CDS <- total_cogs$CDS/100
erba::get_slopePerPhylum(total_cogs, x =  "CDS" , y ="Transcription factors")
erba::get_correlation(total_cogs, "CDS","Transcription factors")

get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(total_cogs, phylum  == x)
  summary(lm(y$`Transcription factors` ~ y$`ORFs(X100)`))
}

sapply(sort(unique(total_cogs$phylum)), get_R2_perPhylum)

####################################################################
# E) COGs Sigma Factors ####

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                 col_names = TRUE,
                                 skip = 2)

archaeas <-  c("Euryarchaeota", "Crenarchaeota")
total_cogs <- total_cogs %>%
  filter(!phylum %in% archaeas)


## plot distribution of transcription factors versus genome size
erba::plot_points(total_cogs, type =  "groups",
                  x = total_cogs$CDS/100,
                  y = total_cogs$`Sigma factors`,
                  filename  =  here::here("figures/F1_cogs_sf_lm_color.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome",
                  xlab = "CDS (x 100)", ymax = 150)

### obtain lm Coefficients per group
total_cogs$CDS <- total_cogs$CDS/100

erba::get_slopePerPhylum(total_cogs, x = "CDS", y = "Sigma factors")
erba::get_correlation(total_cogs, x = "CDS", y = "Sigma factors")

get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(total_cogs, phylum  == x)
  summary(lm(y$`Sigma factors` ~ y$CDS))
}

sapply(sort(unique(total_cogs$phylum)), get_R2_perPhylum)


##########################################################################
# F) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## plot distribution of riboswitches versus genome size

zero_riboswitches <- c("Chlamydiae", "Crenarchaeota")

data_riboswitch <- data_riboswitch %>%
  filter(!phylum %in% zero_riboswitches)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  x = data_riboswitch$CDS/100,
                  y = data_riboswitch$total,
                  filename = here::here("figures/F1_riboswitch_lm_color.tiff"),
                  title = "Transcriptional Riboswitches",
                  xlab = "CDS (x 100)",
                  ylab = "Riboswitches per genome",
                  ymax = 20)

### obtain lm Coefficients per group
data_riboswitch$CDS <- data_riboswitch$CDS/100
erba::get_slopePerPhylum(data_riboswitch, x = "CDS", y = "total")
erba::get_correlation(data_riboswitch,  x = "CDS", y = "total")

get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(data_riboswitch, phylum  == x)
  summary(lm(y$total ~ y$CDS))
}

sapply(sort(unique(data_riboswitch$phylum)), get_R2_perPhylum)


