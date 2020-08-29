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

## plot the number of transcription factors versus sigma factors per genome

TheilSen <- function(..., weights = NULL) {
    mblm::mblm(...)
}

tiff(filename = here::here("figures","F2.tiff"),
     width = 18, height = 12, units = "cm", res = 300)
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
dev.off()

### obtain lm coefficients per phylum
total_cogs %>%
    group_by(phylum) %>%
    summarise(cor = cor(`Sigma factors`, `Transcription factors`))

x <- total_cogs %>% filter(phylum == "Verrucomicrobia")

total_cogs %>%
    group_by(phylum) %>%
    summarise(m = mblm::mblm(sf ~ tf ) )


model <- mblm::mblm(formula = x$`Sigma factors` ~ x$`Transcription factors`)
mblm::mblm(formula = total_cogs$`Sigma factors` ~ total_cogs$`Transcription factors`, dataframe = total_cogs)
