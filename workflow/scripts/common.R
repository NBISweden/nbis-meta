
num_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.numeric))]
    x
}

str_cols <- function(x) {
    x <- x[, unlist(lapply(x, is.character))]
    x
}

process_data <- function(x, output) {
    if (nrow(x)==0) {
        write.table(x, file = output, quote = FALSE, sep="\t")
        quit()
    }
    # Extract only numeric columns
    x_num <- num_cols(x)
    x_num
}
