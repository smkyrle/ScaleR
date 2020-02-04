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
    res <- do.call(cbind, list(RSquare, RMSE)) ## #bind all lists into results
    res <- as.data.frame(res)
    colnames(res) <- c('RSquare', 'RMSE')
  }
  return(res)
}


PLS <- function(x, y){
  NComponents = list() ## define key outputs as lists first to save each for loop iteration
  RSquare <- list() ## ^
  Q2 <- list()
  RMSE <- list()
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
    model <- train( y=training_Y, x= trainingData, method = "pls", ### PLS Model in R using the pls package
                    scale = TRUE, trControl = trainControl("cv", number = 10),
                    tuneLength = 10 )
    predicto <- predict(model, testData)
    RSquare[[paste0("Factor_",i)]] = caret::R2(predicto, testY)
    NComponents[[paste0("Factor_",i)]] = model$bestTune
    res <- do.call(cbind, list(RSquare, NComponents)) ## #bind all lists into results
    res <- as.data.frame(res)
    colnames(res) <- c('RSquare', 'NComponents')
  }
  return(res)
}

