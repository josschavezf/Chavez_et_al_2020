library(erba)

## load data
data_riboswitch_translational <- readr::read_tsv("../Supplementary_Files/Supplementary_Table_15.txt")

## plot distribution of riboswitches versus genome size per phylum
erba::plot_points(data_riboswitch_translational,
                  type = "groups",
                  filename = "figures/translational_riboswitch_lm_color.tiff",
                  title = "Translational Riboswitches",
                  ylab = "Riboswitches per genome",
                  ymax = 30)


# obtain lm Coefficients per group
erba::get_linear_coefficients(data_riboswitch_translational)

# get correlaation index per group
erba::get_correlation(data_riboswitch_translational) %>%
    mutate(cor = round(cor, 2))
