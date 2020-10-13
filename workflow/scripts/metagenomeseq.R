#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(metagenomeSeq)
source("workflow/scripts/common.R")

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

norm <- cbind(str_cols(x), norm)

write.table(x = norm, file = output, quote = FALSE, sep="\t")
