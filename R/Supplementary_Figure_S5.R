library(erba)
library(data.table)
library(dplyr)
library(ggplot2)

# load data

data_cogs <- readxl::read_excel(here::here("data/Table_S3.xlsx"),
                                col_names = TRUE,
                                skip = 2)

data_kos_tf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                            sheet = 1,
                                            col_names = TRUE,
                                            skip = 2)

data_kos_sf <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                  sheet = 2,
                                  col_names = TRUE,
                                  skip = 2)

# select organisms
class_proteobacteria <- c("Alphaproteobacteria",
                          "Betaproteobacteria",
                          "Deltaproteobacteria",
                          "Epsilonproteobacteria",
                          "Gammaproteobacteria")

proteobacteria_cogs <- dplyr::filter(data_cogs, class %in% class_proteobacteria)
proteobacteria_kos_tf <- dplyr::filter(data_kos_tf, class %in% class_proteobacteria)
proteobacteria_kos_sf <- dplyr::filter(data_kos_sf, class %in% class_proteobacteria)

# generate plot function
colors_class <- c("#ed2809", "#ed9209", "#429123", "#3aafd6", "#8515ed", "#ed15d4", "red", "blue", "green")

plot_points_by_class <-
  function(df, x, y, filename = "figure.tiff", plot_title = "", y_title = "", ymax = 200) {
    myplot <- ggplot(df, aes(x , y , colour = class)) +
      geom_point(size = 0.5) +
      ylim(0, ymax) +
      xlim(0,100) +
      geom_smooth(method="lm", se=FALSE, size = 0.5) +
      ggtitle(plot_title) +
      labs(x= "CDS (x 100)",y = y_title) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5,vjust = 0, size = 16, face = "bold"),
            axis.line = element_line(colour = "black", size = 1.2),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.ticks = element_line(size = 1.2, lineend = 2),
            legend.title = element_blank(),
            legend.position = c(0.155,0.85),
            legend.text = element_text(size = 8, hjust = 0),
            legend.key.size = unit(0.45, "cm"),
            legend.spacing.x  = unit(0.05, "cm"),
            legend.margin = unit(0, "cm")
            ) +
      scale_color_manual(values = colors_class, aesthetics = "colour")
    ggsave(filename, plot = myplot, width = 7, height = 5, units = "cm", dpi = 300, scale = 2 )
  }



# plot COG transcriptional regulators for Proteobacteria
plot_points_by_class(proteobacteria_cogs,
                     x = proteobacteria_cogs$CDS/100,
                     y = proteobacteria_cogs$`Transcription factors`,
                     filename = here::here("figures/S5_prot_cogs_tf.tiff"),
                     plot_title = "Transcription Factors",
                     y_title = "Total TFs per genome", ymax = 700)

plot_points_by_class(proteobacteria_cogs,
                     x = proteobacteria_cogs$CDS/100,
                     y = proteobacteria_cogs$`Sigma factors`,
                     filename = here::here("figures/S5_prot_cogs_sf.tiff"),
                     plot_title = "Sigma Factors",
                     y_title = "Total Sigma factors per genome", ymax = 80)


# plot KEGG transcriptional regulators  for Proteobacteria

plot_points_by_class(proteobacteria_kos_tf,
                     x = proteobacteria_kos_tf$CDS/100,
                     y = proteobacteria_kos_tf$total,
                     filename = here::here("figures/S5_prot_kos_tf.tiff"),
                     plot_title = "Transcription Factors",
                     y_title = "Total TFs per genome", ymax = 200)

plot_points_by_class(proteobacteria_kos_sf,
                     x = proteobacteria_kos_sf$CDS/100,
                     y = proteobacteria_kos_sf$total,
                     filename = here::here("figures/S5_prot_kos_sf.tiff"),
                     plot_title = "Sigma Factors",
                     y_title = "Total Sigma factors per genome", ymax = 80)

# Obtain regression model coefficients

# get slope and correlation index

proteobacteria_cogs %>%
  group_by(class) %>%
  mutate(CDS = CDS/100) %>%
  summarise(m = lm(`Transcription factors` ~ CDS)[[1]][2] %>% round(2),
            cor = cor(`Transcription factors`, CDS) %>% round(2))

proteobacteria_cogs %>%
  group_by(class) %>%
  mutate(CDS = CDS/100) %>%
  summarise(m = lm(`Sigma factors` ~ CDS)[[1]][2] %>% round(2),
            cor = cor(`Sigma factors`, CDS) %>% round(2))

proteobacteria_kos_tf %>%
  group_by(class) %>%
  mutate(CDS = CDS/100) %>%
  summarise(m = lm(total ~ CDS)[[1]][2] %>% round(2),
            cor = cor(total, CDS) %>% round(2))

proteobacteria_kos_sf %>%
  group_by(class) %>%
  mutate(CDS = CDS/100) %>%
  summarise(m = lm(total ~ CDS)[[1]][2] %>% round(2),
            cor = cor(total, CDS) %>% round(2))

