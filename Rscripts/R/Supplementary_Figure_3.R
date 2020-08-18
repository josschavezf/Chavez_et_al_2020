library(dplyr)
library(erba)
library(here)
library(readxl)
library(reshape2)

##################################################################

# KEGG

## load data
total_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                 col_names = TRUE,
                                 sheet = 1,
                                 skip = 2)
total_kos_sigma <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                   col_names = TRUE,
                                   sheet = 2,
                                   skip = 2)
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## summarise data
total_kos_tf <- total_kos_tf %>%
  select("organism","phylum", "ORFs(X100)","total")
colnames(total_kos_tf)[4] <- "Transcription factors"

total_kos_sigma <- total_kos_sigma %>%
  select("organism","total")
colnames(total_kos_sigma)[2] <- "Sigma factors"

data_riboswitch <- data_riboswitch %>%
  select("organism","total")
colnames(data_riboswitch)[2] <- "Riboswitches"

data_merged <- merge(total_kos_tf, total_kos_sigma)
data_merged <- merge(data_merged, data_riboswitch)

data_melted <- reshape2::melt(data_merged,
                              id.vars = c("organism","phylum","ORFs(X100)"),
                              variable.name = "type",
                              value.name = "total")

## plot data
sapply(unique(data_melted$phylum), function(x) {data_melted %>%
    filter(phylum == x) %>%
    erba::plot_regulators(paste0("S3", x),x) } )
