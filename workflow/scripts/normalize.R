#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

num_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.numeric))]
    x
}

str_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.character))]
    x
}

method <- snakemake@params$method

# Read the counts
x <- read.delim(snakemake@input[[1]], row.names = 1)
# Remove unclassified
x <- x[row.names(x)!="Unclassified", ]
# Extract only numeric columns
x_num <- num_cols(x)

if (method %in% c("TMM", "RLE")) {
    library(edgeR, quietly = TRUE)
    # Create DGE
    obj <- DGEList(x_num)
    # Calculate norm factors
    obj <- calcNormFactors(obj, method = method)
    # Calculate cpms
    norm <- cpm(obj, normalized.lib.sizes = TRUE)
} else {
    library(metagenomeSeq, quietly = TRUE)
    obj <- newMRexperiment(x_num)
    norm <- MRcounts(obj, norm = TRUE)
}

norm <- cbind(str_cols(x), norm)

write.table(x = norm, file = snakemake@output[[1]], quote = FALSE, sep="\t")
