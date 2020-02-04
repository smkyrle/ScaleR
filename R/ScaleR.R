#' @title  ScaleR: An R package for investigating scaling factors effects on models
#'
#' @description This package employs Ridge regression or PLS over multiple scaling factors.
#'
#' @param x is a data matrix of predictors, features in columns and samples in rows.
#'
#' @parma y is the outcome variable.
#'
#' @param inter is the Scalaing factor intervals, e.g 0.1, 0.05, 0.01.
#'
#' @param methods is the model type you want to employ (methods='Ridge), default is PLS.
#'
#' @return matrix of correspongin model assesment metrics
#'
#' @examples ScaleR(x,y,inter=0.1, meethod='PLS')
#'
#' @export











ScaleR <- function(x,y, method, inter ,plots){
  Scaling_Factor = seq(0,1, inter)
  Data_Sets <- SCALE(x)

  if (method=='Ridge') {
    print('Employing Ridge Regression')
    Ridge(Data_Sets, y)}
  else {
    print('Employing PLS')
    res<- PLS(Data_Sets, y) }
  return(res)}
