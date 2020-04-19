log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Libload
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Read in homer file
homer.df <- readr::read_tsv(snakemake@input[[1]])
# Remove the feature details
homer.df <- homer.df %>% mutate(Annotation = sapply(Annotation, function(x){strsplit(as.character(x), split = " (", fixed = TRUE)[[1]][1]}))
# Keep only required columns
homer.df <- homer.df %>% select(c("Chr", "Start", "End", "Annotation", "Gene Name", "Gene Type", "Distance to TSS", "Gene Alias", "Gene Description"))
# Write output
write.table(homer.df, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE, na = "")

sessionInfo()
