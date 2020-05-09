# attach erba
library(erba)

##################################################################

# KEGG

# load data
total_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")
total_kos_sigma <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_9.txt")
data_riboswitch <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_10.txt")

# create sumarized table

data_summarized <- total_kos_tf %>%
  select(organism, ORFs, phylum, class) %>%
  bind_cols(Transcription_Factor = total_kos_tf$total) %>%
  bind_cols(Sigma_Factor = total_kos_sigma$total) %>%
  bind_cols(Riboswitches = data_riboswitch$total)

data_summarized <- reshape2::melt(data_summarized, id.vars = c("organism","phylum", "class","ORFs"))

colnames(data_summarized)[5] <- "regulator"
colnames(data_summarized)[6] <- "total"

# create function to plot ####
colors_reg <- c("red", "blue", "darkgreen")
names(colors_reg) <- c("Sigma_Factor", "Transcription_Factor", "Riboswitches")

plot_regulators <- function(data, filename = "figure.tiff", title = "") {
  tiff(filename, width = 1234, height = 880, units = 'px', res = 100)
  myplot <- ggplot(data, aes(x = ORFs, y = total, colour = regulator)) +
    geom_point() +
    #ylim(0,ymax) +
    #xlim(0,100) +
    geom_smooth(method="lm", se=FALSE) +
    ggtitle(title) +
    labs(x= "ORFs (x 100)",y = "Transcriptional regulators per genome") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
          axis.line = element_line(colour = "black", size = 1.5),
          axis.title = element_text(size = 26, face = "bold"),
          axis.text = element_text(size = 22, face = "bold"),
          axis.ticks = element_line(size = 1.5, lineend = 2),
          legend.title = element_blank(),
          legend.position = c(0.165,0.9),
          legend.text = element_text(size = 20, face = "bold"),
          legend.key.size = unit(0.45, "in")) +
    scale_color_manual(values = colors_reg,  aesthetics = "colour")
  print(myplot)
  dev.off()
}

# plot regulators vs ORFs per group ####
for (i in selected_phylogeny) {
  filename = paste0("figures/sup2_",i, ".tiff")
  x <- data_summarized %>% filter(phylum == i)
  plot_regulators(x, title = i, filename)
}
