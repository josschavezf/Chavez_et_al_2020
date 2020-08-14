# phylum-specific Transcription factors

library(erba)

# import data
data_kos_tf <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_7.txt")

# perform fisher test
fisher_tf <- erba::enrichment_fisher(data_kos_tf,
                                     column_names = kos_tf)

# import KO families description
kos_tf_names <- readr::read_tsv("../Supplementary_files/Supplementary_Table_2.txt")

# replace column names
colnames(fisher_tf[[2]]) <- tf_names$Prefix

# set data structure
library(reshape2)
tf_melted <- as_tibble(melt(fisher_tf[[2]]))

colnames(tf_melted) <- c("phylum", "KO", "value")

# filter phylum-specific families
tf_filtered <- tf_melted %>%
  filter(value == "Inf") %>%
  arrange(phylum) %>%
  select(phylum, KO)

# Add family description and write table
inner_join(tf_filtered, kos_tf_names) %>%
  readr::write_tsv("../Supplementary_files/Supplementary_Table_13.txt")
