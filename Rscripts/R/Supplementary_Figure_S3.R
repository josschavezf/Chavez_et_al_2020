library(dplyr)
library(here)
library(readxl)

##########################################################################
# KEGG

## load data
transcription_factors <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                            sheet = 1,
                                            col_names = TRUE,
                                            skip = 2)

sigma_factors <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                    sheet = 2,
                                    col_names = TRUE,
                                    skip = 2)

kos_regulators <- transcription_factors[1:4]
kos_regulators$tf <- transcription_factors$total
kos_regulators$sf <- sigma_factors$total

## plot the number of transcription factors versus sigma factors per genome

TheilSen <- function(..., weights = NULL) {
    mblm::mblm(...)
}

tiff(filename = here::here("figures","S2.tiff"),
     width = 18, height = 12, units = "cm", res = 300)
ggplot(kos_regulators, aes(x = tf, y = sf)) +
    geom_point(size = 0.2) +
    geom_smooth(method = "TheilSen", se = FALSE, size = 0.5) +
    labs(x = "Transcription factors",
         y = "Sigma factors") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 8, face = "bold")) +
    facet_wrap(~phylum, scales = "free")
dev.off()
