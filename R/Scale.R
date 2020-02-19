SCALE <- function(d, Scaling_Factor) {
  Data_Sets <- list()
  for (i in Scaling_Factor){
    x.centered <- apply(d, 2, function(x) x - mean(x)) ### x-mean(x)
    # Then we perform scaling on the mean-centered matrix
    x.sc <- apply(x.centered, 2, function(x) x/(sd(x)^i)) ## new x/sd raised to the power of sclaing factor
    Data_Sets[[paste0("Factor_", i)]] <- x.sc
  }
  return(Data_Sets)
}
