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
lm(total_cogs$`Transcription Factors` ~ total_cogs$`ORFs(X100)`)
cor(total_cogs$`Transcription Factors`, total_cogs$`ORFs(X100)`) %>% round(2)

## plot quadratic distribution of transcription factors versus genome size
erba::plot_points(total_cogs, type = "general",
                  column_orfs = total_cogs$`ORFs(X100)`,
                  column_total = total_cogs$`Transcription factors`,
                  filename = here::here("figures/cogs_tf_factor_quadratic.tiff"),
                  title = "Transcription factors",
                  ylab = "Transcription factors per genome", ymax = 800,
                  model.degree = 2)

y <- total_cogs$`Transcription factors`
x <- total_cogs$`ORFs(X100)`
x2 <- x^2

summary(lm(y ~ x))
summary(lm(y ~ x  + x2))


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
cor(total_cogs$`Sigma factors`, total_cogs$`ORFs(X100)`) %>% round(2)


## plot quadratic distribution of sigma factors versus genome size
erba::plot_points(total_cogs, type =  "general",
                  column_orfs = total_cogs$`ORFs(X100)`,
                  column_total = total_cogs$`Sigma factors`,
                  filename =  here::here("figures/cogs_sigma_factor_quadratic.tiff"),
                  title ="Sigma factors",
                  ylab = "Sigma factors per genome", ymax = 150,
                  model.degree = 2)

y <- total_cogs$`Sigma factors`
x <- total_cogs$`ORFs(X100)`
x2 <- x^2

summary(lm(y ~ x))
summary(lm(y ~ x + x2))
cor(y, x)

##########################################################################
# C) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)
colnames(data_riboswitch)[2] <- "ORFs"

## plot distribution of riboswitches versus genome size
erba::plot_exception(data_riboswitch,
                     filename = here::here("figures/riboswitch_lm.tiff"),
                     title = "Transcriptional Riboswitches",
                     ylab = "Riboswitches per genome",
                     ymax = 80,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)

data_riboswitch  %>%
    filter(phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)

data_riboswitch  %>%
    filter(!phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)


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
erba::get_correlation(total_cogs, x = "ORFs(X100)", y = "Transcription Factors")
erba::get_slopePerPhylum(total_cogs, x =  "ORFs(X100)" , y ="Transcription Factors")


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

##########################################################################
# F) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## plot distribution of riboswitches versus genome size
zero_groups <- data_riboswitch %>%
  group_by(phylum) %>%
  summarise(total = sum(total)) %>%
  filter(total == 0) %>%
  select(phylum) %>% unlist()

data_riboswitch <- data_riboswitch %>%
  filter(!phylum %in% zero_groups)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  column_total = data_riboswitch$total,
                  column_orfs = data_riboswitch$`ORFs(X100)`,
                  filename = "figures/riboswitch_lm_color.tiff",
                  title = "Transcriptional Riboswitches",
                  ylab = "Riboswitches per genome",
                  ymax = 80)

### obtain lm Coefficients per group
erba::get_correlation(data_riboswitch, x = "ORFs(X100)", y = "total")
erba::get_slopePerPhylum(data_riboswitch, x = "ORFs(X100)", y = "total")
