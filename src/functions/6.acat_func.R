# Using function and getting a list
acat_column <- function(x) {
  # Convert non-NA values to numeric
  Pvals <- as.numeric(x[!is.na(x)])
  
  # Call ACAT() with numeric Pvals
  o <- ACAT(Pvals = Pvals)
  return(o)
}
output <- lapply(output_table, function(col) {
  if (sum(!is.na(col)) > 0) {
    acat_column(col)
  } else {
    NA
  }
})

# Directly usng ACAT to get a dataframe
output_df <- do.call(rbind, lapply(output_table, function(col) {
  if (sum(!is.na(col)) > 0) {
    pvalue <- as.numeric(na.omit(col))
    o <- ACAT(Pvals = pvalue)
    return(data.frame(ACAT_result = o))
  } else {
    return(data.frame(ACAT_result = NA))
  }
}))

