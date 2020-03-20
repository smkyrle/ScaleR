#' @title  ScaleR: An R package for investigating scaling factors effects on models
#'
#' @description This package employs Ridge regression or PLS over multiple scaling factors. Default is PLS
#'
#' @param x is a data matrix of predictors, features in columns and samples in rows.
#'
#' @parma y is the outcome variable.
#'
#' @param inter is the Scalaing factor intervals, e.g 0.1, 0.05, 0.01. For instance
#'
#' @param methods is the model type you want to employ (methods='Ridge), default is PLS.
#'
#' @param plot this determines if you rewuire a plot of R" and Q2 against scaling factor, default is TRUE
#'
#' @return matrix of correspongin model assesment metrics
#'
#' @examples ScaleR(x,y,inter=0.1, meethod='PLS')
#'
#' @export




ScaleR <- function(x,y, method='PLS', inter=NULL, plot=TRUE){

  ### Test train 80/20 split
  if (is.null(inter)){print('interval not provided, using default value of 0.1')
    inter=0.1}
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  if (nrow(x)==0){
    print('Must provide x as a numeric Data Set')
  }
  if (nrow(x)!=nrow(y)){
    print('Must provide x and y with same number of samples/rows')
  }
  y<- as.matrix(y)
  set.seed(123)
  trainingIndex <- sample(1:nrow(x), 0.8*nrow(x)) # indices for 80% training data ### can edit this/ take the training and test spilit out of for loop incase we need directly comparable results
  trainingData <- x[trainingIndex, ] # training data
  testData <- x[-trainingIndex, ] # test data
  trainingY <- y[trainingIndex, ] ## note the use of index for the split nrow xtrain= nrow ytrain ##
  testY <- y[-trainingIndex, ] # test data
  inter <- inter
  Scaling_Factor <- seq(0,1, as.numeric(inter))### Scaling_Factor
  Data_Sets <- SCALE(trainingData, Scaling_Factor) ### Data_Sets of initial training with diffrent scalings

  if (method=='Ridge') {
    print('Employing Ridge Regression')
    Res <- Ridge(Data_Sets, trainingY)
    Res$Scaling_Factor <- Scaling_Factor
    if (plot){
      plot(as.numeric(Res$Scaling_Factor),as.numeric(Res$RSquared_Y),type="l",col="red", xlim=c(0, 1), ylim=c(0, max(as.numeric(Res$RSquared_Y)+0.1)), xlab='', ylab='')
      par(new=TRUE)
      lines(as.numeric(Res$Scaling_Factor),as.numeric(Res$QSquared_Y), col='green')
      title(main="Scaling Factor vs R2Y and Q2Y",
            xlab="Scaling Factor ", ylab="R2Y & Q2Y")
      par(xpd=TRUE)
      legend("topright", title="Metric",
             c("R2","Q2"), fill=c('red','green'))}} ###

  if (method=='PLS') {
    print('Employing PLS')
    Res <- PLS(Data_Sets, trainingY)
    Res$Scaling_Factor <- Scaling_Factor
    if (plot){
      plot(as.numeric(Res$Scaling_Factor),as.numeric(Res$RSquared_Y),type="l",col="red", xlim=c(0, 1), ylim=c(0, max(as.numeric(Res$RSquared_Y)+0.1)), xlab='', ylab='')
      par(new=TRUE)
      lines(as.numeric(Res$Scaling_Factor),as.numeric(Res$QSquared_Y), col='green')
      title(main="Scaling Factor vs R2Y and Q2Y",
            xlab="Scaling Factor ", ylab="R2Y & Q2Y")
      par(xpd=TRUE)
      legend("topright", title="Metric",
             c("R2","Q2"), fill=c('red','green'))}} ###
  Res <- list(Res, testY, trainingY )
  names(Res) <- c('Results', 'Test_Y', 'Training_Y')
  return(Res)
}
