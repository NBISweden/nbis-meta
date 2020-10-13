
num_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.numeric))]
    x
}

str_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.character))]
    x
}

process_data <- function(input, output) {
    # Read the counts
    x <- read.delim(input, row.names = 1, sep = "\t", header = TRUE)
    # Remove unclassified
    if ("Unclassified" %in% row.names(x)){
        x <- x[row.names(x)!="Unclassified", ]
    }
    if (nrow(x)==0) {
        write.table(x, file = output, quote = FALSE, sep="\t")
        quit()
    }

    # Extract only numeric columns
    x_num <- num_cols(x)
    x_num
}
