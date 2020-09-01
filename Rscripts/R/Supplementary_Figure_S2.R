library(dplyr)
library(erba)
library(here)
library(readxl)

##########################################################################
# A) Riboswitch #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 3,
                                      col_names = TRUE,
                                      skip = 2)

colnames(data_riboswitch)[2] <- "ORFs"

zero_riboswitches <- c("Chlamydiae", "Crenarchaeota")
data_riboswitch <- data_riboswitch %>% filter(phylum != zero_riboswitches)

## plot distribution of riboswitches versus genome size
erba::plot_exception(data_riboswitch,
                     filename = here::here("figures/S2_riboswitch_copies_lm.tiff"),
                     title = "Transcriptional Riboswitches",
                     ylab = "Riboswitches per genome",
                     ymax = 80,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)

data_riboswitch  %>%
    filter(phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)

data_riboswitch  %>%
    filter(!phylum == "Firmicutes")  %>%
    summarise(cor = cor(total, ORFs)) %>%
    round(2)

##########################################################################
# B) Riboswitch by phylum #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 3,
                                      col_names = TRUE,
                                      skip = 2)

## plot distribution of riboswitches versus genome size
zero_groups <- data_riboswitch %>%
  group_by(phylum) %>%
  summarise(total = sum(total)) %>%
  filter(total == 0) %>%
  select(phylum) %>% unlist()

data_riboswitch <- data_riboswitch %>% filter(!phylum %in% zero_groups)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  column_total = data_riboswitch$total,
                  column_orfs = data_riboswitch$`ORFs(X100)`,
                  filename = "figures/S2_riboswitch_copies_lm_color.tiff",
                  title = "Transcriptional Riboswitches",
                  ylab = "Riboswitches per genome",
                  ymax = 80)

### obtain lm Coefficients per group
erba::get_slopePerPhylum(data_riboswitch, x = "ORFs(X100)", y = "total")
erba::get_correlation(data_riboswitch, x = "ORFs(X100)", y = "total")
