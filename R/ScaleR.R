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
#' @param seed_val this will allow user to define a seed number if not provided will use the default 1234 
#'
#'@ param k is the number of best perforing models according to Robustness of Cross Validation (R2Y/Q2Y)
#'
#' @return matrix of correspongin model assesment metrics
#'
#' @examples ScaleR(x,y,inter=0.1, meethod='PLS')
#'
#' @export

ScaleR <- function(x,y, method='PLS', inter=NULL, plot=TRUE, seed_val=1234, k=NULL){
    
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
    if (is.null(seed_val)){
        seed_val=1234}
    seed_val <- as.numeric(seed_val)
    str(seed_val)
    y<- as.matrix(y)
    trainingData <- as.matrix(x)
    inter <- inter     
    Scaling_Factor <- seq(0,1, as.numeric(inter))### Scaling_Factor 
    Data_Sets <- SCALE(trainingData, Scaling_Factor)
    method_opt= 'RCV' ### Data_Sets of initial training with diffrent scalings 
    
    

    if (method=='Ridge') {
                print('Employing Ridge Regression')
                Res <- Ridge(Data_Sets, y, seed_val)
                Res$Scaling_Factor <- Scaling_Factor
                Res$RCV <- as.numeric(Res$RSquared_Y)/as.numeric(Res$QSquared_Y)
                B.ind <- which.minn(unlist(Res[method_opt]), k=2)
                Res <- as.list(Res)
                if (as.numeric(Scaling_Factor[B.ind[1]])==0)  {
                    inter_2 = as.numeric(Scaling_Factor[B.ind[2]])/10 } 
                if (as.numeric(Scaling_Factor[B.ind[1]])==1){
                    inter_2 = (as.numeric(Scaling_Factor[B.ind[1]]) - as.numeric(Scaling_Factor[B.ind[2]]))/10} else {
                    inter_2 = as.numeric(Scaling_Factor[B.ind[1]])/10}
                if (as.numeric(Scaling_Factor[B.ind[1]]) < as.numeric(Scaling_Factor[B.ind[2]])) {
                    Scaling_Factor_B <- seq(as.numeric(Scaling_Factor[B.ind[1]]),as.numeric(Scaling_Factor[B.ind[2]]), as.numeric(inter_2))}else{ 
                    Scaling_Factor_B <- seq(as.numeric(Scaling_Factor[B.ind[2]]),as.numeric(Scaling_Factor[B.ind[1]]), as.numeric(inter_2))}
                Data_Sets2 <- SCALE(trainingData, Scaling_Factor_B)
                Res2 <- Ridge(Data_Sets2, y, seed_val)
               
                Res2$Scaling_Factor <- Scaling_Factor_B
                Res2$RCV <- as.numeric(Res2$RSquared_Y)/as.numeric(Res2$QSquared_Y)
                Res3 <- mapply(c,Res[1:6], as.list(Res2[1:6]), SIMPLIFY=FALSE)
                list.2d <- sapply(Res3[1:6], cbind)
                df.2d <- as.data.frame(list.2d)
                #df.2d <- list.2d[order(-as.numeric(Res3$RSquared_Y)),]
                df.2d <- df.2d[!duplicated(df.2d$Scaling_Factor),]
                df.2d <- df.2d[order(-as.numeric(df.2d$Scaling_Factor)),]
                Res[1:6] <- df.2d
               
                
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
        Data.k <- SCALE(trainingData, as.numeric(Res$Scaling_Factor[k.RCV]))
        model <- PLS.k(Data.k, y, seed_val)
        Res$model <- model}} ### 
    
            

    if (method=='PLS') {
                print('Employing PLS')
                Res <- PLS(Data_Sets, y, seed_val)
                Res$Scaling_Factor <- Scaling_Factor
                Res$RCV <- as.numeric(Res$RSquared_Y)/as.numeric(Res$QSquared_Y)
                B.ind <- which.minn(unlist(Res[method_opt]), k=2)
                if is.null(B.ind)){print('All RCV(s) are negative, using R-squared')
                method_opt=RSquared_Y
                B.ind <- which.minn(unlist(Res[method_opt]), k=2)} 
                Res <- as.list(Res)
                if (as.numeric(Scaling_Factor[B.ind[1]])==0)  {
                    inter_2 = as.numeric(Scaling_Factor[B.ind[2]])/10 } 
                if (as.numeric(Scaling_Factor[B.ind[1]])==1){
                    inter_2 = (as.numeric(Scaling_Factor[B.ind[1]]) - as.numeric(Scaling_Factor[B.ind[2]]))/10} else {
                    inter_2 = as.numeric(Scaling_Factor[B.ind[1]])/10}
                if (as.numeric(Scaling_Factor[B.ind[1]]) < as.numeric(Scaling_Factor[B.ind[2]])) {
                    Scaling_Factor_B <- seq(as.numeric(Scaling_Factor[B.ind[1]]),as.numeric(Scaling_Factor[B.ind[2]]), as.numeric(inter_2))}else{ 
                    Scaling_Factor_B <- seq(as.numeric(Scaling_Factor[B.ind[2]]),as.numeric(Scaling_Factor[B.ind[1]]), as.numeric(inter_2))}
                Data_Sets2 <- SCALE(trainingData, Scaling_Factor_B)
                Res2 <- PLS(Data_Sets2, y, seed_val)
               
                Res2$Scaling_Factor <- Scaling_Factor_B
                Res2$RCV <- as.numeric(Res2$RSquared_Y)/as.numeric(Res2$QSquared_Y)
                Res3 <- mapply(c,Res[1:6], as.list(Res2[1:6]), SIMPLIFY=FALSE)
                list.2d <- sapply(Res3[1:6], cbind)
                df.2d <- as.data.frame(list.2d)
                #df.2d <- list.2d[order(-as.numeric(Res3$RSquared_Y)),]
                df.2d <- df.2d[!duplicated(df.2d$Scaling_Factor),]
                df.2d <- df.2d[order(-as.numeric(df.2d$Scaling_Factor)),]
                Res[1:6] <- df.2d
               

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
            Data.k <- SCALE(trainingData, as.numeric(Res$Scaling_Factor[k.RCV]))
            model <- PLS.k(Data.k, y, seed_val)
            Res$model <- model}} ### 
return(Res)
}
