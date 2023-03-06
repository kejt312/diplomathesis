#! /usr/bin/Rscript

library(Maaslin2)

features <- file.path('/','tmp','features.tsv')
metadata <- file.path('/','tmp','metadata.tsv')

fit_data <- Maaslin2(
    features, metadata, file.path('/','tmp','maaslin2'), 
    fixed_effects = c('diagnosis'),
    random_effects = c('study'),
    normalization = 'NONE',
    transform = 'NONE'
)

