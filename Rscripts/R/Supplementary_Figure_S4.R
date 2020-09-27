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
      geom_point(size = 0.05) +
      ylim(0, ymax) +
      xlim(0,100) +
      geom_smooth(method="lm", se=FALSE, size = 0.2) +
      ggtitle(plot_title) +
      labs(x= "ORFs (x 100)",y = y_title) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5,vjust = 0, size = 8, face = "bold"),
            axis.line = element_line(colour = "black", size = 0.4),
            axis.title = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 6, face = "bold"),
            axis.ticks = element_line(size = 0.2, lineend = 2),
            legend.title = element_blank(),
            legend.position = c(0.155,0.85),
            legend.text = element_text(size = 5.5, hjust = 0),
            legend.key.size = unit(0.25, "cm"),
            legend.spacing.x  = unit(0.05, "cm"),
            legend.margin = unit(0, "cm")
            ) +
      scale_color_manual(values = colors_class, aesthetics = "colour")
    ggsave(filename, plot = myplot, width = 9, height = 6.5, units = "cm", dpi = 300 )
  }



# plot COG transcriptional regulators for Proteobacteria
plot_points_by_class(proteobacteria_cogs,
                     x = proteobacteria_cogs$`ORFs(X100)`,
                     y = proteobacteria_cogs$`Transcription factors`,
                     filename = here::here("figures/S4_prot_cogs_tf.tiff"),
                     plot_title = "Transcription Factors",
                     y_title = "Total TFs per genome", ymax = 700)

plot_points_by_class(proteobacteria_cogs,
                     x = proteobacteria_cogs$`ORFs(X100)`,
                     y = proteobacteria_cogs$`Sigma factors`,
                     filename = here::here("figures/S4_prot_cogs_sf.tiff"),
                     plot_title = "Sigma Factors",
                     y_title = "Total Sigma factors per genome", ymax = 80)


# plot KEGG transcriptional regulators  for Proteobacteria

plot_points_by_class(proteobacteria_kos_tf,
                     x = proteobacteria_kos_tf$`ORFs(X100)`,
                     y = proteobacteria_kos_tf$total,
                     filename = here::here("figures/S4_prot_kos_tf.tiff"),
                     plot_title = "Transcription Factors",
                     y_title = "Total TFs per genome", ymax = 200)

plot_points_by_class(proteobacteria_kos_sf,
                     x = proteobacteria_kos_sf$`ORFs(X100)`,
                     y = proteobacteria_kos_sf$total,
                     filename = here::here("figures/S4_prot_kos_sf.tiff"),
                     plot_title = "Sigma Factors",
                     y_title = "Total Sigma factors per genome", ymax = 80)

# Obtain regression model coefficients

# get slopes

proteobacteria_cogs %>%
  group_by(class) %>%
  summarise(m = lm(`Transcription factors` ~ `ORFs(X100)`)[[1]][2] %>% round(2) )

proteobacteria_cogs %>%
  group_by(class) %>%
  summarise(m = lm(`Sigma factors` ~ `ORFs(X100)`)[[1]][2] %>% round(2) )

proteobacteria_kos_tf %>%
  group_by(class) %>%
  summarise(m = lm(total ~ `ORFs(X100)`)[[1]][2] %>% round(2))

proteobacteria_kos_sf %>%
  group_by(class) %>%
  summarise(m = lm(total ~ `ORFs(X100)`)[[1]][2] %>% round(2))


# Obtain correlation index
get_correlation_class <- function(x) {
  classes <- x %>%
    distinct(class) %>%
    arrange(class) %>%
    unlist()
  for (i in classes) {
    print(i)
    table <- filter(x, class == i)
    cor(table$total, table$ORFs) %>%
      round(2) %>%
      print()
  }
}

get_correlation_class(proteobacteria_cogs_tf)
get_correlation_class(proteobacteria_cogs_sigma)
get_correlation_class(proteobacteria_kos_tf)
get_correlation_class(proteobacteria_kos_sigma)

