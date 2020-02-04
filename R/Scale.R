SCALE <- function(x) {
  Data_Sets <- list()
  for (i in Scaling_Factor){
    x.centered <- apply(x, 2, function(x) x - mean(x))
    # Then we perform scaling on the mean-centered matrix
    x.sc <- apply(x.centered, 2, function(x) x/(sd(x)^i))
    Data_Sets[[paste0("Factor_", i)]] <- x.sc
  }
  return(Data_Sets)
}
