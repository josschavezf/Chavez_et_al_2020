library(erba)
library(data.table)
library(dplyr)
library(ggplot2)

# load data
data_cogs_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_6.txt")
data_cogs_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_8.txt")

data_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")
data_kos_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_9.txt")
data_riboswitch <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_10.txt")

# select organisms
class_proteobacteria <- c("Alphaproteobacteria",
                          "Betaproteobacteria",
                          "Deltaproteobacteria",
                          "Epsilonproteobacteria",
                          "Gammaproteobacteria",
                          "Zetaproteobacteria")

proteobacteria_cogs_tf <- filter(data_cogs_tf, class %in% class_proteobacteria)
proteobacteria_cogs_sigma <- filter(data_cogs_sigma, class %in% class_proteobacteria)
proteobacteria_kos_tf <- filter(data_kos_tf, class %in% class_proteobacteria)
proteobacteria_kos_sigma <- filter(data_kos_sigma, class %in% class_proteobacteria)
proteobacteria_riboswitch <- filter(data_riboswitch, class %in% class_proteobacteria)

# generate plot function
colors_class <- c("#ed2809", "#ed9209", "#429123", "#3aafd6", "#8515ed", "#ed15d4", "red", "blue", "green")
plot_points_by_class <-
  function(x, tiff_file = "figure.tiff", plot_title = "", y_title = "", ymin = 0, ymax = 200) {
  tiff(tiff_file, width = 1234, height = 880, units = 'px', res = 100)
  myplot <- ggplot(x, aes(x = ORFs, y = total, colour = class)) +
    geom_point() +
    ylim(ymin, ymax) +
    xlim(0,100) +
    geom_smooth(method="lm", se=FALSE) +
    ggtitle(plot_title) +
    labs(x= "ORFs (x 100)",y = y_title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
          axis.line = element_line(colour = "black", size = 1.5),
          axis.title = element_text(size = 26, face = "bold"),
          axis.text = element_text(size = 22, face = "bold"),
          axis.ticks = element_line(size = 1.5, lineend = 2),
          legend.title = element_blank(),
          legend.position = c(0.155,0.82),
          legend.text = element_text(size = 20),
          legend.key.size = unit(0.45, "in") ) +
    scale_color_manual(values = colors_class, aesthetics = "colour")
  print(myplot)
  dev.off()
}

# filter exception within Deltaproteobacteria
deltaproteobacteria_riboswitch_filtered <- proteobacteria_riboswitch %>%
  filter(class == "Deltaproteobacteria") %>%
  filter(total > 1)
proteobacteria_riboswitch_filtered <- proteobacteria_riboswitch %>%
  filter(!class == "Deltaproteobacteria") %>%
  bind_rows(deltaproteobacteria_riboswitch_filtered)

# plot transcriptional regulators by COGs for Proteobacteria
plot_points_by_class(proteobacteria_cogs_tf, "figures/cogs_proteobacteria_tf.tiff", "Transcription Factors", "Total TFs per genome", ymax = 700)
plot_points_by_class(proteobacteria_cogs_sigma, "figures/cogs_proteobacteria_sigma.tiff", "Sigma Factors", "Total sigma factors per genome", ymax = 80)

# plot transcriptional regulators by KOs for Proteobacteria
plot_points_by_class(proteobacteria_kos_tf, "figures/kos_proteobacteria_tf.tiff", "Transcription Factors", "Total TFs per genome")
plot_points_by_class(proteobacteria_kos_sigma, "figures/kos_proteobacteria_sigma.tiff", "Sigma Factors", "Total sigma factors per genome", ymax = 80)
plot_points_by_class(proteobacteria_riboswitch_filtered, "figures/proteobacteria_riboswitch.tiff", "Riboswitches", "Total riboswitches per genome", ymax = 30)

# Obtain regression model coefficients

# define function

get_slope_by_class <- function(x) {
  classes <- x %>% distinct(class) %>% unlist()
  for (i in classes) {
    print(i)
    filter(x, class == i) %>%
      lm(formula = total~ORFs) %>%
      print()
  }
}

# obtain slopes
get_slope_by_class(proteobacteria_cogs_tf)
get_slope_by_class(proteobacteria_cogs_sigma)

get_slope_by_class(proteobacteria_kos_tf)
get_slope_by_class(proteobacteria_kos_sigma)
get_slope_by_class(proteobacteria_riboswitch_filtered)


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
get_correlation_class(proteobacteria_riboswitch_filtered)

