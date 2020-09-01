library(dplyr)
library(here)
library(readxl)
library(ggplot2)

##########################################################################
# COG

## load data
total_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                 col_names = TRUE,
                                 skip = 2)

colnames(total_cogs)[6:7] <- c("tf", "sf")
colnames(total_cogs)[5] <- "ORFs"

## remove archaeas

archaeas <- c("Crenarchaeota", "Euryarchaeota")
total_cogs <- dplyr::filter(total_cogs, !phylum %in% archaeas)

## reorder by phylum
phylum_order <- sort(unique(total_cogs$phylum))[-c(3,8)]
phylum_order[8:9] <- c("Chlamydiae" , "Tenericutes" )

total_cogs$phylum <- factor(total_cogs$phylum, ordered = TRUE, levels = phylum_order)

## plot the number of transcription factors versus sigma factors per genome

TheilSen <- function(..., weights = NULL) {
    mblm::mblm(...)
}

ggplot(total_cogs, aes(x = tf, y = sf)) +
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
ggsave(filename = here::here("figures","F2.tiff"), device = "tiff",
       width = 9, height = 6, units = "cm", dpi = 300, scale = 2)

### obtain lm coefficients per phylum

slopeThielsen <- function(x) {
    df <- total_cogs %>% dplyr::filter(phylum == x)
    m <- mblm::mblm(sf ~ tf ,df)[[1]][2] %>% round(2)
    names(m) <- x
    return(m)
}

sapply(levels(total_cogs$phylum), slopeThielsen)

model <- mblm::mblm(formula = x$`Sigma factors` ~ x$`Transcription factors`)
mblm::mblm(formula = total_cogs$`Sigma factors` ~ total_cogs$`Transcription factors`, dataframe = total_cogs)
