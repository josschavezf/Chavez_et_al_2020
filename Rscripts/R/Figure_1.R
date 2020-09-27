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
                  column_orfs = total_cogs$`ORFs(X100)`,
                  column_total = total_cogs$`Transcription factors`,
                  filename = here::here("figures/cogs_tf_factor_lm.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome", ymax = 800)

### obtain lm Coefficients in general
lm(total_cogs$`Transcription factors` ~ total_cogs$`ORFs(X100)`)
cor(total_cogs$`Transcription factors`, total_cogs$`ORFs(X100)`)

summary(lm(total_cogs$`Transcription factors` ~ total_cogs$`ORFs(X100)`))

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
                  column_orfs = total_cogs$`ORFs(X100)`,
                  column_total = total_cogs$`Sigma factors`,
                  filename =  here::here("figures/cogs_sigma_factor_lm.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 150)

## obtain lm Coefficients in general
lm(total_cogs$`Sigma factors` ~ total_cogs$`ORFs(X100)` )
cor(total_cogs$`Sigma factors` , total_cogs$`ORFs(X100)` )

summary(lm(total_cogs$`Sigma factors` ~ total_cogs$`ORFs(X100)`))
##########################################################################
# C) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

colnames(data_riboswitch)[2] <- "ORFs"

zero_riboswitches <- c("Chlamydiae", "Crenarchaeota")
data_riboswitch <- data_riboswitch %>% dplyr::filter(!phylum %in% zero_riboswitches)

erba::plot_exception(data_riboswitch,
                     filename = here::here("figures/F1_riboswitch_lm.tiff"),
                     title = "Transcriptional Riboswitches",
                     ylab = "Riboswitches per genome",
                     ymax = 15,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)

filter(data_riboswitch, phylum == "Firmicutes") %>%
  summarise(cor = cor(total, ORFs) )
filter(data_riboswitch, !phylum == "Firmicutes") %>%
  summarise(cor = cor(total, ORFs) )

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
                  column_total = total_cogs$`Transcription Factors`,
                  column_orfs = total_cogs$`ORFs(X100)`,
                  filename = here::here("figures/cogs_tf_factor_lm_color.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome", ymax = 800)

## obtain lm Coefficients per group
erba::get_correlation(total_cogs, "ORFs(X100)","Transcription factors")
erba::get_slopePerPhylum(total_cogs, x =  "ORFs(X100)" , y ="Transcription Factors")

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
                  column_total = total_cogs$`Sigma factors`,
                  column_orfs = total_cogs$`ORFs(X100)`,
                  filename  =  here::here("figures/cogs_sigma_factor_lm_color.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 150)

### obtain lm Coefficients per group
erba::get_correlation(total_cogs, x = "ORFs(X100)", y = "Sigma factors")
erba::get_slopePerPhylum(total_cogs, x = "ORFs(X100)", y = "Sigma factors")

get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(total_cogs, phylum  == x)
  summary(lm(y$`Sigma factors` ~ y$`ORFs(X100)`))
}

sapply(sort(unique(total_cogs$phylum)), get_R2_perPhylum)


##########################################################################
# F) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)
colnames(data_riboswitch)[2] <- "ORFs"

## plot distribution of riboswitches versus genome size

zero_riboswitches <- c("Chlamydiae", "Crenarchaeota")

data_riboswitch <- data_riboswitch %>%
  filter(!phylum %in% zero_riboswitches)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  x = data_riboswitch$ORFs,
                  y = data_riboswitch$total,
                  filename = here::here("figures/F1_riboswitch_lm_color.tiff"),
                  title = "Transcriptional Riboswitches",
                  xlab = "ORFs (x 100)",
                  ylab = "Riboswitches per genome",
                  ymax = 15)

### obtain lm Coefficients per group
erba::get_slopePerPhylum(data_riboswitch, x = "ORFs", y = "total")
erba::get_correlation(data_riboswitch,  x = "ORFs", y = "total")
get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(data_riboswitch, phylum  == x)
  summary(lm(y$total ~ y$ORFs))
}

sapply(sort(unique(data_riboswitch$phylum)), get_R2_perPhylum)


