library(erba)
library(circlize)
library(ComplexHeatmap)

# set vector with colors
col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# plot_legend ####
lgd = Legend(col_fun = col_fun, title = "log10(odd.ratio)", legend_height = unit(8, "cm"))

tiff(filename = "figures/heatmap_legend.tiff",height = 1200, res = 300)
draw(lgd)
dev.off()

###############################################################################
# COGs Transcription Factors ####

data_cogs_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_11.txt")

## enrichment analysis ####
fisher_cogs_tf <- enrichment_fisher(data_cogs_tf,
                                       column_names = cogs_tf,
                                       phylogeny = selected_phylogeny)

## replace colnames ####
enrichment_cogs_tf <- log10(fisher_cogs_tf[[2]] + 1e-5)
cogs_tf_names <- readr::read_tsv("../Supplementary_files/Supplementary_Table_1.txt")
colnames(enrichment_cogs_tf) <- as.character(cogs_tf_names$Prefix)

## plot heatmap ####
tiff(filename = "figures/heatmap_cogs_tf.tiff",width = 800,  units = 'px', res = 100)
myplot <- Heatmap(enrichment_cogs_tf,
                  name = "log10(odd.ratio)",
                  col = col_fun,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 14),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 14),
                  show_column_names = FALSE,
                  show_heatmap_legend = FALSE,
                  height = unit(10, "cm"))
print(myplot)
dev.off()

###############################################################################
# COGs Sigma Factors ####

data_cogs_sigma <- readxl::read_excel(here::here("data/Table_S6.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## enrichment analysis ####
archaeas <- c("Crenarchaeota","Euryarchaeota")
fisher_cogs_sigma <- erba::enrichment_fisher(data_cogs_sigma,
                                       column_names = cogs_sigma,
                                       phylogeny = selected_phylogeny[!selected_phylogeny %in% archaeas])
enrichment_cogs_sigma <- log10(fisher_cogs_sigma[[2]] + 1e-5)

## replace colnames ####
cogs_sigma_names <- readxl::read_excel(here::here("data/Table_S1.xlsx"),
                                       sheet = 3,
                                       col_names = TRUE,
                                       skip = 2)
cogs_sigma_names <- cogs_sigma_names %>%
    select(COG, Prefix) %>%
    distinct()
colnames(enrichment_cogs_sigma) <- as.character(cogs_sigma_names$Prefix)

## plot heatmap ####
tiff(filename = here::here("figures/F3_heatmap_cogs_sf.tiff"), height = 800, units = 'px', res = 100)
myplot <- Heatmap(enrichment_cogs_sigma,  name = "log10(odd.ratio)", col = col_fun,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 14),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 14),
                  show_heatmap_legend = FALSE,
                  height = unit(10, "cm"))
print(myplot)
dev.off()


###############################################################################
# KOs Transcription Factors ####

data_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")

## enrichment analysis ####
fisher_kos_tf <- enrichment_fisher(data_kos_tf,
                                    column_names = kos_tf,
                                    phylogeny = selected_phylogeny,
                                   replace_by = 3000)

## replace colnames ####
enrichment_kos_tf <- log10(fisher_kos_tf[[2]] + 1e-5)
kos_tf_names <- readr::read_tsv("../Supplementary_files/Supplementary_Table_2.txt")

colnames(enrichment_kos_tf) <- as.character(kos_tf_names$Prefix)

## plot heatmap ####
tiff(filename = "figures/heatmap_kos_tf.tiff", width = 800,  units = 'px', res = 100)
myplot <- Heatmap(enrichment_kos_tf,  name = "log10(odd.ratio)", col = col_fun,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 14),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 14),
                  show_column_names = FALSE,
                  show_heatmap_legend = FALSE,
                  height = unit(10, "cm"))
print(myplot)
dev.off()

###############################################################################

# KOs Sigma Factors ####

data_kos_sigma <- readxl::read_excel(here::here("data/Table_S4.xlsx"),
                                     sheet = 2,
                                     col_names = TRUE,
                                     skip = 2)

## enrichment analysis ####
fisher_kos_sigma <- enrichment_fisher(data_kos_sigma,
                                      column_names = kos_sigma,
                                      phylogeny = selected_phylogeny[!selected_phylogeny %in% archaeas],
                                      replace_by = 1000)
enrichment_kos_sigma <- log10(fisher_kos_sigma[[2]] + 1e-5)

## replace colnames ####
kos_sigma_names <- readxl::read_excel(here::here("data/Table_S1.xlsx"),
                                      sheet = 4,
                                      col_names = TRUE,
                                      skip = 2)
colnames(enrichment_kos_sigma) <- as.character(kos_sigma_names$Prefix)

## plot heatmap ####
tiff(filename = here::here("figures/F3_heatmap_kos_sf.tiff"), height = 600, width = 600, units = 'px', res = 100)
myplot <- Heatmap(enrichment_kos_sigma,  name = "log10(odd.ratio)", col = col_fun,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 14),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 14),
                  show_heatmap_legend = FALSE,
                  height = unit(10, "cm"),
                  width = unit(10, "cm"))
print(myplot)
dev.off()


###############################################################################
# unique TF kos for archaea

# import data
data_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")

# calculate total regulators per phylum per KO
data_kos_tf_summarised <-  data_kos_tf %>%
  group_by(phylum) %>%
  summarise_if(is.numeric, sum)

# transpose counts to keep kos as rows and phylum as columns
data_kos_tf_summarised <-  tblhelpr::transpose_tibble(data_kos_tf_summarised,
                                                      col_names = "phylum",
                                                      id_col = "KO")

# filter data where archaeas have elements
data_archaea <- filter(data_kos_tf_summarised,
                       Crenarchaeota > 0 | Euryarchaeota > 0 ) %>%
  mutate(total = rowSums(.[-1])) %>%
  mutate(total_archaea = Crenarchaeota + Euryarchaeota)

# filter data from archaea to eliminate kos present in other phyla
data_archaea_unique <- filter(data_archaea,
                              total == total_archaea)

# import kos description
kos_tf_names <- readr::read_tsv("../Supplementary_files/Supplementary_Table_2.txt")

# see description for kos unique for archaea
filter(kos_tf_names, KO %in% data_archaea_unique$KO)


###############################################################################

# Riboswitch #####

data_riboswitch <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

## enrichment analysis ####
id_transcriptional <- colnames(data_riboswitch)[grep("RF[0-9]",colnames(data_riboswitch))]

fisher_riboswitch <- erba::enrichment_fisher(data_riboswitch,
                                             column_names = id_transcriptional,
                                             replace_by = 1000,
                                             phylogeny = erba::selected_phylogeny)
riboswitch <- log10(fisher_riboswitch[[2]] + 1E-5)

## replace colnames ####
riboswitch_names <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                       sheet = 1,
                                       col_names = TRUE,
                                       skip = 2)
riboswitch_names <- dplyr::filter(riboswitch_names, Rfam %in% id_transcriptional)
colnames(riboswitch) <- riboswitch_names$Description

## plot heatmap ####
tiff(filename = "figures/F3_heatmap_riboswitch.tiff", width = 1234, height = 880, units = 'px', res = 100)
myplot <- Heatmap(riboswitch,  name = "log10(odd.ratio)", col = col_fun,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 16),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 16),
                  show_heatmap_legend = FALSE)
print(myplot)
dev.off()
