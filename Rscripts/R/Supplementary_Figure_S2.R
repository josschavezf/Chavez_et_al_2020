library(dplyr)
library(erba)
library(here)
library(readxl)
library(reshape2)

##################################################################

# COG

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                 col_names = TRUE,
                                 skip = 2)
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## summarise data
data_merged <- merge(total_cogs, data_riboswitch[c("organism", "total")])
colnames(data_merged)[8] <- c("Riboswitches")

data_melted <- reshape2::melt(data_merged,
                              id.vars = c("organism","phylum", "class","Specie","ORFs(X100)"),
                              variable.name = "type",
                              value.name = "total")

## plot data
sapply(unique(data_melted$phylum), function(x) {data_melted %>%
         filter(phylum == x) %>%
         erba::plot_regulators(paste0("S2", x),x) } )
