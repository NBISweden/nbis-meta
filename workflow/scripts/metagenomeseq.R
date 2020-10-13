#!/usr/bin/env Rscript

library(metagenomeSeq)
source("common.R")

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

method <- snakemake@params$method
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]

x_num <- process_data(input, output)

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
