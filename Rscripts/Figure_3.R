library(erba)

## load data ####
attach(erba::tf_repressor_non_repressor)

## plot transcription factors as repressor or non repressor function ####
plot_repressor(tf_repressor_non_repressor,
               filename = "figures/tf_repressor_non_repressor.tiff",
               title = "Transcription factors",
               ylab = "Transcription factors per genome")

### obtain correlation index for repressors and non repressors ####
coef(lm(non_repressor~ORFs, data = tf_repressor_non_repressor))
coef(lm(repressor~ORFs, data = tf_repressor_non_repressor))

tf_repressor_non_repressor  %>%
    summarise(cor = cor(non_repressor, ORFs)) %>%
    round(2)

tf_repressor_non_repressor  %>%
    summarise(cor = cor(repressor, ORFs)) %>%
    round(2)
