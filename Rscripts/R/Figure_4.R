###############################################################################
# comparative dotplot
library(tidyverse)
library(reshape2)

# read riboswitch Rfam names and description
riboswitch_names <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                       sheet = 1,
                                       col_names = TRUE,
                                       skip = 2)

# load riboswitch count data
data_riboswitch_transcriptional <- readxl::read_excel(here::here("data/Table_S5.xlsx"),
                                      sheet = 2,
                                      col_names = TRUE,
                                      skip = 2)

data_riboswitch_translational <- readxl::read_excel(here::here("data/Table_S9.xlsx"),
                                                      col_names = TRUE,
                                                      skip = 2)

# filter riboswitch names present within each data
id_transcriptional <- colnames(data_riboswitch_transcriptional)[grep("RF[0-9]",colnames(data_riboswitch_transcriptional))]

id_translational <- colnames(data_riboswitch_translational)[grep("RF[0-9]",colnames(data_riboswitch_translational))]

riboswitch_names_transcriptional <- dplyr::filter(riboswitch_names, Rfam %in% id_transcriptional)

riboswitch_names_translational <- dplyr::filter(riboswitch_names, Rfam %in% id_translational)

# replace Rfam names by riboswitch names
colnames(data_riboswitch_transcriptional)[5:36] <- riboswitch_names_transcriptional$Description
colnames(data_riboswitch_translational)[5:39] <- riboswitch_names_translational$Description

# melt tables to get riboswitch names as rows
info.vars <- c("organism", "ORFs(X100)", "phylum", "class")

df_riboswitch_transcriptional <- reshape2::melt(data_riboswitch_transcriptional[-37],
                                      id.vars = info.vars,
                                      variable.name = "riboswitch",
                                      value.name = "transcriptional")

df_riboswitch_translational <- reshape2::melt(data_riboswitch_translational[-40],
                                    id.vars = info.vars,
                                    variable.name = "riboswitch",
                                    value.name = "translational")

# merge transcriptional and translational tables

df_riboswitch <- full_join(df_riboswitch_transcriptional,
                           df_riboswitch_translational)

df_riboswitch[is.na(df_riboswitch)] <- 0

# count the number of organisms per phylum
count_organisms <- data_riboswitch_transcriptional %>%
  group_by(phylum) %>%
  count()

# generate the summarized table
df_riboswitch_group <- df_riboswitch %>%
  group_by(phylum, riboswitch) %>%
  summarise(transcriptional = sum(transcriptional),
            translational = sum(translational)) %>%
  mutate(total = translational+transcriptional) %>%
  mutate(transcriptional_percent = transcriptional*100/total) %>%
  right_join(count_organisms) %>%
  mutate(mean_transcriptional = transcriptional/n) %>%
  mutate(mean_translational = translational/n) %>%
  mutate(mean_riboswitches_per_organism = (transcriptional+translational)/n) %>%
  mutate(mean_transcriptional_percent = mean_transcriptional*100/mean_riboswitches_per_organism)

# replace total == 0 by NA
df_riboswitch_group$total <- na_if(df_riboswitch_group$total, 0)
df_riboswitch_group$mean_riboswitches_per_organism <- na_if(df_riboswitch_group$mean_riboswitches_per_organism,0)

## reorder by phylum
phylum_order <- c("Firmicutes",
                  "Tenericutes",
                  "Actinobacteria",
                  "Proteobacteria",
                  "Spirochaetes",
                  "Bacteroidetes",
                  "Planctomycetes",
                  "Verrucomicrobia",
                  "Chlamydiae",
                  "Euryarchaeota",
                  "Crenarchaeota")

df_riboswitch_group$phylum <- factor(df_riboswitch_group$phylum,
                                     ordered = TRUE, levels = phylum_order)

ribos_order <- c("TPP", levels(df_riboswitch_group$riboswitch)[-2] )

df_riboswitch_group$riboswitch <- factor(df_riboswitch_group$riboswitch,
                                     ordered = TRUE, levels = ribos_order)


# plot points as mean percent of transcriptional riboswitches

ggplot(df_riboswitch_group,
       aes(x = riboswitch, y = phylum,
           size = mean_riboswitches_per_organism, color = mean_transcriptional_percent)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 6, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 6),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "#f5f5f5"),
        legend.position = "bottom",
        legend.margin = margin( r = 0.2, unit = "cm"),
        legend.title = element_text(size = 6, hjust = 0.5),
        legend.text = element_text(size = 6, hjust = 0.5),
        legend.key = element_blank(),
        legend.spacing.x = unit(0.05, 'cm')) +
  guides(size = guide_legend(title = "Mean of riboswitches per organism",
                             title.position = "top",
                             nrow = 1,
                             order = 1),
         colour = guide_colorbar(title = "% of transcriptional riboswitches",
                                 title.position = "top",
                                 order = 2,
                                 barwidth = 7)
  ) +
  scale_color_gradientn(colours = rev(rainbow(5))) +
  scale_y_discrete(limits = rev(levels(as.factor(df_riboswitch_group$phylum))) ) +
  scale_size_continuous(range = c(0.5,5), breaks = seq(0.1,1, 0.1) )
ggsave(filename = here::here("figures","F42.tiff"), device = "tiff",
       width = 8, height = 5, units = "cm", dpi = 300, scale = 2)

