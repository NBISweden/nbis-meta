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
if ("Unclassified" %in% row.names(x)){
    x <- x[row.names(x)!="Unclassified", ]
}
# Extract only numeric columns
x_num <- num_cols(x)

if (method %in% c("TMM", "RLE")) {
    library(edgeR)
    # Create DGE
    obj <- DGEList(x_num)
    # Calculate norm factors
    obj <- calcNormFactors(obj, method = method)
    # Calculate cpms
    norm <- cpm(obj, normalized.lib.sizes = TRUE)
} else if (method == "RPKM") {
    library(edgeR)
    # Extract gene length column
    gene_length <- x_num$Length
    x_num <- x_num[, colnames(x_num) != "Length"]
    obj <- DGEList(x_num)
    # Calculate norm factors
    obj <- calcNormFactors(obj, method = "TMM")
    # Calculate RPKM
    norm <- rpkm(obj, normalized.lib.sizes = TRUE, gene.length = gene_length)
} else {
    library(metagenomeSeq)
    obj <- newMRexperiment(x_num)
    norm <- MRcounts(obj, norm = TRUE)
}

norm <- cbind(str_cols(x), norm)

write.table(x = norm, file = snakemake@output[[1]], quote = FALSE, sep="\t")
