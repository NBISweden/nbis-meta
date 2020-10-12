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
x <- read.delim(snakemake@input[[1]], row.names = 1, sep = "\t", header = TRUE)
# Remove unclassified
if ("Unclassified" %in% row.names(x)){
    x <- x[row.names(x)!="Unclassified", ]
}
if (nrow(x)==0) {
    write.table(x, file = snakemake@output[[1]], quote = FALSE, sep="\t")
    quit()
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
    smat = lapply(1:ncol(x_num), function(i) {
	    sort(x_num[which(x_num[, i]>0),i], decreasing = TRUE)
	})
    if(any(sapply(smat,length)==1)) {
        fh <-file(snakemake@output[[1]])
        writeLines(c("WARNING: Sample with one or zero features",
                     "Cumulative Sum Scaling failed for sample"), fh)
        close(fh)
        quit()
    }
    norm <- MRcounts(obj, norm = TRUE)
}

norm <- cbind(str_cols(x), norm)

write.table(x = norm, file = snakemake@output[[1]], quote = FALSE, sep="\t")
