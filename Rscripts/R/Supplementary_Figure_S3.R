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

## remove archaeas

archaeas <- c("Crenarchaeota", "Euryarchaeota")
kos_regulators <- dplyr::filter(kos_regulators, !phylum %in% archaeas)

## reorder by phylum
phylum_order <- sort(unique(kos_regulators$phylum))[-c(3,8)]
phylum_order[8:9] <- c("Chlamydiae" , "Tenericutes" )

kos_regulators$phylum <- factor(kos_regulators$phylum, ordered = TRUE, levels = phylum_order)


## plot the number of transcription factors versus sigma factors per genome

TheilSen <- function(..., weights = NULL) {
    mblm::mblm(...)
}


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
ggsave(filename = here::here("figures","S3.tiff"), device = "tiff",
       width = 9, height = 6, units = "cm", dpi = 300, scale = 2)

### obtain lm coefficients per phylum

slopeThielsen <- function(x) {
    df <- kos_regulators %>% dplyr::filter(phylum == x)
    m <- mblm::mblm(sf ~ tf ,df)[[1]][2] %>% round(2)
    names(m) <- x
    return(m)
}

sapply(levels(kos_regulators$phylum), slopeThielsen)
