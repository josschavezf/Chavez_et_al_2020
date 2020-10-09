
library(dplyr)
library(VennDiagram)

# Supplementary Figure 1
##########################################################################

# read data
df_cog <- read.table("table_cog.txt", sep = "\t")
colnames(df_cog) <- c("gene","org","COG", "n1", "n2")

df_ko <- read.table("table_ko.txt", sep = "\t")
colnames(df_ko) <- c("gene","org", "KO")



# read cog and ko ids for transcripcion factors and sigma factors
cogs_tf <- readxl::read_excel("data/Table_S1.xlsx",
                              sheet = 1,
                              skip = 2)
kos_tf <- readxl::read_excel("data/Table_S1.xlsx",
                             sheet = 2,
                             skip = 2)

cogs_sf <- readxl::read_excel("data/Table_S1.xlsx",
                              sheet = 3,
                              skip = 2)

kos_sf <- readxl::read_excel("data/Table_S1.xlsx",
                              sheet = 4,
                              skip = 2)

# filter data

df_cog_tf <- df_cog %>% filter(COG %in% cogs_tf$COG)
df_ko_tf <- df_ko %>% filter(KO %in% kos_tf$KO)

df_cog_sf <- df_cog %>% filter(COG %in% unique(cogs_sf$COG) )
df_ko_sf <- df_ko %>% filter(KO %in% kos_sf$KO)

# Plot Venn diagrams

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(x = list(COG = unique(df_cog_tf$gene), KO = unique(df_ko_tf$gene)),
             filename = "figures/S1_tf.tiff",
             fill = myCol[1:2],
             fontface = "bold",
             cat.fontface = "bold",
             main.fontface = "bold",
             cex = 0.7,
             cat.cex = 0.7,
             height = 6,
             width = 6,
             units = "cm",
             resolution = 300)

venn.diagram(x = list(COG = unique(df_cog_sf$gene), KO = unique(df_ko_sf$gene)),
             filename = "figures/S1_sf.tiff",
             fill = myCol[1:2],
             fontface = "bold",
             cat.fontface = "bold",
             main.fontface = "bold",
             cex = 0.6,
             cat.cex = 0.7,
             height = 6,
             width = 6,
             units = "cm",
             resolution = 300)

