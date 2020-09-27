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
data_riboswitch <- data_riboswitch %>% filter(!phylum %in% zero_riboswitches)

## plot distribution of riboswitches versus genome size
erba::plot_exception(data_riboswitch,
                     filename = here::here("figures/S2_riboswitch_copies_lm.tiff"),
                     title = "Transcriptional Riboswitches",
                     ylab = "Copies of riboswitches per genome",
                     ymax = 80,
                     exception_group = "Firmicutes")

### obtain lm Coefficients in general
lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)
lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs)[[1]][2] %>% round(2)

filter(data_riboswitch, phylum == "Firmicutes") %>%
  summarise(cor(total, ORFs) )
filter(data_riboswitch, !phylum == "Firmicutes") %>%
  summarise(cor(total, ORFs) )

summary(lm(data = data_riboswitch[data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs) )
summary(lm(data = data_riboswitch[!data_riboswitch$phylum =="Firmicutes",],
   formula = total~ORFs) )

##########################################################################
# B) Riboswitch by phylum #####

## load data
data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 3,
                                      col_names = TRUE,
                                      skip = 2)
colnames(data_riboswitch)[2] <- "ORFs"

## plot distribution of riboswitches versus genome size
zero_groups <- data_riboswitch %>%
  group_by(phylum) %>%
  summarise(total = sum(total)) %>%
  filter(total == 0) %>%
  select(phylum) %>% unlist()

data_riboswitch <- data_riboswitch %>% filter(!phylum %in% zero_groups)

erba::plot_points(data_riboswitch,
                  type = "groups",
                  x = data_riboswitch$ORFs,
                  y = data_riboswitch$total,
                  filename = "figures/S2_riboswitch_copies_lm_color.tiff",
                  title = "Transcriptional Riboswitches by phyla",
                  xlab = "ORFs (x100)",
                  ylab = "Copies of riboswitches per genome",
                  ymax = 80)

### obtain lm Coefficients per group
erba::get_slopePerPhylum(data_riboswitch, x = "ORFs", y = "total")
erba::get_correlation(data_riboswitch, x = "ORFs", y = "total")


get_R2_perPhylum <- function(x) {
  y <- dplyr::filter(data_riboswitch, phylum  == x)
  summary(lm(y$total ~ y$ORFs))
}

sapply(sort(unique(data_riboswitch$phylum)), get_R2_perPhylum)
