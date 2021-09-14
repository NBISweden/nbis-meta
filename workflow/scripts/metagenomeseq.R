#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(metagenomeSeq)
source("workflow/scripts/common.R")

method <- snakemake@params$method
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
normalize <- TRUE
# Read the counts
x <- read.delim(input, row.names = 1, sep = "\t", header = TRUE)
# Get sample names
sample_names <- colnames(x)[unlist(lapply(x, is.numeric))]
# Extract row names
xrownames <- row.names(x)[row.names(x)!="Unclassified"]
# Get info names
info_names <- colnames(x)[unlist(lapply(x, is.character))]

# Remove unclassified
if ("Unclassified" %in% row.names(x)){
    x <- x[row.names(x)!="Unclassified", ]
}

x <- as.data.frame(x, row.names=xrownames)
colnames(x) <- append(info_names, sample_names)

# Returns a vector in the case of 1 sample only
x_num <- process_data(x, output)

# If only one sample, set normalize=FALSE
if (length(sample_names) == 1) {
    print("ONLY ONE SAMPLE! WILL NOT RUN CSS")
    normalize <- FALSE
    x_num <- as.data.frame(x_num, row.names = rownames(x))
    colnames(x_num) <- sample_names
}

# Turn data into new experimentobject
obj <- newMRexperiment(x_num)
too_few_features <- FALSE
# In cases with only one sample, check whether there are enough features to run
# CSS
if (is.null(ncol(x_num))) {
    if (length(x_num) <= 1) {
        too_few_features <- TRUE
    }
} else { # Do the same type of check for many samples
    smat <- lapply(1:ncol(x_num), function(i) {
        sort(x_num[which(x_num[, i]>0),i], decreasing = TRUE)
    })
    if (any(sapply(smat,length)==1)) {
        too_few_features <- TRUE
    }
}
if (too_few_features == TRUE) {
    fh <-file(snakemake@output[[1]])
    writeLines(c("WARNING: Sample with one or zero features",
                 "Cumulative Sum Scaling failed for sample"), fh)
    close(fh)
    quit()
}

# Normalize
norm <- MRcounts(obj, norm = normalize)

# Add info columns back
norm <- cbind(str_cols(x), norm)
colnames(norm) <- append(info_names, sample_names)

# Convert to numeric
norm <- as.data.frame(norm)

# Set sample names
colnames(norm)[unlist(lapply(norm, is.numeric))] <- sample_names

# Write output
write.table(x = norm, file = output, quote = FALSE, sep="\t")
