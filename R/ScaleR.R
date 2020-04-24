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
#' @param plot this determines if you require a plot of R2Y and Q2Y against scaling factor, default is TRUE
#'
#' @param seed this will allow user to define a seed number if not provided will use the default 1234 
#'
#'@ param k is the number of best perforing models according to Robustness of Cross Validation (R2Y/Q2Y)
#'
#' @return matrix of correspongin model assesment metrics
#'
#' @examples ScaleR(x,y,inter=0.1, meethod='PLS')
#'
#' @export


ScaleR <- function(x,y, method='PLS', inter=NULL, plot=TRUE, seed=NULL, k=NULL){
    
    ## Initital checks inter, data sets etc. ## 
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
    ## set seed to default if null/not privided##
    if (is.null(seed)){
        seed=1234}
    seed <- as.numeric(seed)
    y<- as.matrix(y)
    trainingData <- as.matrix(x)
    inter <- inter     
    Scaling_Factor <- seq(0,1, as.numeric(inter))### Scaling_Factor 
    Data_Sets <- SCALE(trainingData, Scaling_Factor) ### Data_Sets of initial training with diffrent scalings 
    str(Data_Sets)
    

    if (method=='Ridge') {
                print('Employing Ridge Regression')
                Res <- Ridge(Data_Sets, y)
                Res$Scaling_Factor <- Scaling_Factor
                Res$RCV <- as.numeric(Res$RSquared_Y)/as.numeric(Res$QSquared_Y)
                
    if (plot){
            plot(as.numeric(Res$Scaling_Factor),as.numeric(Res$RSquared_Y),type="l",col="red", xlim=c(0, 1), ylim=c(0, max(as.numeric(Res$RSquared_Y)+0.1)), xlab='', ylab='')
            par(new=TRUE)
            lines(as.numeric(Res$Scaling_Factor),as.numeric(Res$QSquared_Y), col='green')
            title(main="Scaling Factor vs R2 and Q2",
            xlab="Scaling Factor ", ylab="R2 & Q2")
            par(xpd=TRUE)
            legend("topright", title="Metric",
            c("R2","Q2"), fill=c('red','green'))}### 
    
    if (k>0 & method=='Ridge'){
            k.RCV <- Ind.RCV.k(Res$RCV, k)
            Data.k <- Data_Sets[c(k.RCV)]
            str(Data.k)
            model <- Ridge.k(Data.k, y)
            Res <- as.list(Res)
            Res$model <- model}}
    
            

    if (method=='PLS') {
                print('Employing PLS')
                Res <- PLS(Data_Sets, y)
                Res$Scaling_Factor <- Scaling_Factor
                Res$RCV <- as.numeric(Res$RSquared_Y)/as.numeric(Res$QSquared_Y)

    if (plot){
            plot(as.numeric(Res$Scaling_Factor),as.numeric(Res$RSquared_Y),type="l",col="red", xlim=c(0, 1), ylim=c(0, max(as.numeric(Res$RSquared_Y)+0.1)), xlab='', ylab='')
            par(new=TRUE)
            lines(as.numeric(Res$Scaling_Factor),as.numeric(Res$QSquared_Y), col='green')
            title(main="Scaling Factor vs R2 and Q2",
            xlab="Scaling Factor ", ylab="R2 & Q2")
            par(xpd=TRUE)
            legend("topright", title="Metric",
            c("R2","Q2"), fill=c('red','green'))}
        
        if (k>0 & method=='PLS'){
            k.RCV <- Ind.RCV.k(Res$RCV, k)
            Data.k <- Data_Sets[c(k.RCV)]
            str(Data.k)
            model <- PLS.k(Data.k, y)
            Res <- as.list(Res)
            Res$model <- model}} ### 
return(Res)
}
