library(erba)

attach(erba::tf_repressor_non_repressor)

# define function
library(ggplot2)
plot_repressor_by_phylum <- function(data, filename, title) {
  tiff(filename = filename, width = 1234, height = 880, units = 'px', res = 100)
  myplot <- ggplot(data) +
    geom_point(aes(x = ORFs, y = repressor, colour = "repressor")) +
    geom_point(aes(x = ORFs, y = non_repressor, colour = "non_repressor")) +
    geom_abline(aes(intercept = intercept(ORFs, repressor, slope(ORFs,repressor)), slope = slope(ORFs,repressor), color = "repressor")) +
    geom_abline(aes(intercept = intercept(ORFs, non_repressor, slope(ORFs,non_repressor)), slope = slope(ORFs,non_repressor), color = "non_repressor")) +
    ggtitle(title) +
    labs(x= "ORFs (x 100)", y = "Transcription factors per genome") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
          axis.title = element_text(size = 26, face="bold"),
          axis.text = element_text(size = 22, face="bold"),
          axis.line = element_line(colour = "black", size=1.5),
          axis.ticks = element_line(size = 1.5, lineend = 2),
          legend.title = element_blank(),
          legend.position = c(0.13,0.85),
          legend.text = element_text(size = 20),
          legend.key.size = unit(0.45, "in")) +
    scale_color_manual(values = colors_represor_activator,  aesthetics = "colour")
  print(myplot)
  dev.off()

}

# plot repressor and non repressor by phlyum
for (i in selected_phylogeny) {
  df <- filter(tf_repressor_non_repressor, phylum == i)
  plot_repressor_by_phylum(df,
                           filename = paste0("figures/repressor_",i,".tiff"),
                           title = i)
}


# get linear coeffients by phylum for non_repressor
library(dplyr)
tf_repressor_non_repressor %>%
  group_by(phylum) %>%
  summarise(slope = lm(non_repressor~ORFs)$coefficients[2])

# get correlation index
tf_repressor_non_repressor %>%
  group_by(phylum) %>%
  summarise(cor = round(cor(non_repressor, ORFs), 2) )

# get linear coeffients by phylum for repressor
tf_repressor_non_repressor %>%
  group_by(phylum) %>%
  summarise(slope = lm(repressor~ORFs)$coefficients[2])

# get correlation index
tf_repressor_non_repressor %>%
  group_by(phylum) %>%
  summarise(cor = round(cor(repressor, ORFs), 2) )
