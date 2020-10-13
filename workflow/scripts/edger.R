#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
source("workflow/scripts/common.R")

library(edgeR)
method <- snakemake@params$method
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]

# Read the counts
x <- read.delim(input, row.names = 1, sep = "\t", header = TRUE)
# Remove unclassified
if ("Unclassified" %in% row.names(x)){
    x <- x[row.names(x)!="Unclassified", ]
}

x_num <- process_data(x, output)

if (method %in% c("TMM", "RLE")) {
    # Create DGE
    obj <- DGEList(x_num)
    # Calculate norm factors
    obj <- calcNormFactors(obj, method = method)
    # Calculate cpms
    norm <- cpm(obj, normalized.lib.sizes = TRUE)
} else if (method == "RPKM") {
    # Extract gene length column
    gene_length <- x_num$Length
    x_num <- x_num[, colnames(x_num) != "Length"]
    obj <- DGEList(x_num)
    # Calculate norm factors
    obj <- calcNormFactors(obj, method = "TMM")
    # Calculate RPKM
    norm <- rpkm(obj, normalized.lib.sizes = TRUE, gene.length = gene_length)
}

norm <- cbind(str_cols(x), norm)

write.table(x = norm, file = output, quote = FALSE, sep="\t")
