Ridge <-function(x, y){
  RSquare <- list()
  # Q2 <- list()
  RMSE <- list()
  lambda <- list()
  res <- list()
  for (i in seq_along(x)) {
    set.seed(123)
    inputData <- x[i]
    inputData <- as.data.frame(inputData)
    trainingIndex <- sample(1:nrow(inputData), 0.8*nrow(inputData)) # indices for 80% training data ### can edit this/ take the training and test spilit out of for loop incase we need directly comparable results
    trainingData <- inputData[trainingIndex, ] # training data
    testData <- inputData[-trainingIndex, ] # test data
    training_Y <- y[trainingIndex, ] ## note the use of index for the split nrow xtrain= nrow ytrain ##
    testY <- y[-trainingIndex, ]
    lambdas <- 10^seq(3, -2, by = -.1) ### glmnet has some bult in lambda range. Here we have values ranging from λ=1010 to λ=10−2, essentially covering the full range of scenarios from the null model containing only the intercept, to the least squares fit.
    cv.out = cv.glmnet(x=data.matrix(trainingData), y=training_Y, alpha = 0) # Fit ridge regression model on training data, default 10 fold crossvalidation
    bestlam = cv.out$lambda.min  # Select lamda that minimizes training MSE
    predictions <- predict(cv.out, s = bestlam, data.matrix(testData))  # predict on test data
    RootMSE = caret::RMSE(predictions, testY)  ## RMSE test
    Rsquared = caret::R2(predictions, testY)   ## R-Squared test
    RSquare[[paste0("Factor_", i)]] <- Rsquared ## need to figure out what is going on with R2 and Q2 in R
    RMSE[[paste0("Factor_", i)]] <- RootMSE
    lambda[[paste0("Factor_", i)]] <- cv.out$lambda.min
    res <- do.call(cbind, list(RSquare, RMSE, lambda)) ## #bind all lists into results
    res <- as.data.frame(res)
    colnames(res) <- c('RSquare', 'RMSE', 'lambda')
  }
  return(res)
}


PLS <- function(x, y){
  NComponents = list() ## define key outputs as lists first to save each for loop iteration
  RSquare <- list() ## ^
  QSquare <- list()
  model_out <- list()
  res <- list()
  CV_RMSEP <- list()
  for (i in seq_along(x)) {
    set.seed(123)
    train.data <- x[[i]] ### One of the training Data sets extrpolated
    model <- pls::plsr(as.matrix(trainingY)~as.matrix(train.data), data=as.data.frame(train.data), validation = "CV") ## validation-cv, 10 fold cross val
    CVRMSEP <- pls::RMSEP(model,'CV')$val ### cross-model-validated RMSEP calculated from the training data set- as not to leak info
    NComp <- which.min(CVRMSEP) ## index of minimun CVRMSEP to detemine optimal number of components
    RSquare[[paste0(names(x[i]))]] <- pls::R2(model, 'train')$val[as.numeric(NComp)] ## Take R2 for optimal no. comps
    QSquare[[paste0(names(x[i]))]]  <- pls::R2(model, 'CV')$val[as.numeric(NComp)] ## Q2 also known as cross-validated R2
    CV_RMSEP[[paste0("Factor_",i)]]<- min(CVRMSEP)
    NComponents[[paste0("Factor_",i)]]<- as.numeric(NComp)## save no.componenst
    model_out[[paste0("Factor_",i)]] <- model ### save each model iteritviely
    res <- do.call(cbind, list(RSquare, QSquare, NComponents, CV_RMSEP, model_out)) ## #bind all lists into results
    res <- as.data.frame(res)
    colnames(res) <- c('RSquared', 'QSquared', 'NComponents', 'CV_RMSEP', 'model')
  }
  return(res)
}

