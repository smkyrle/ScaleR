#### Ridge Function
Ridge <-function(x, y){
  RSquare <- list()
  CV_RMSE <- list()
  lambda_out <- list()
  res <- list()
  QSquare <- list()
  for (i in seq_along(x)) {
    set.seed(123)
    train.data <- as.matrix(x[[i]]) ### One of the training Data sets extrpolated
    lambdas <- 10^seq(3, -2, by = -.1) ### glmnet has some bult in lambda range. Here we have values ranging from λ=1010 to λ=10−2, essentially covering the full range of scenarios from the null model containing only the intercept, to the least squares fit.
    cv.out = glmnet::cv.glmnet(x=train.data, y=as.matrix(y), lambda=lambdas, keep=TRUE, type.measure="mse", alpha = 0, standardize= FALSE) # Fit ridge regression model on training data, default 10 fold crossvalidation
    bestlam = cv.out$lambda.min  # Select lamda that minimizes training MSE
    mse <- min(cv.out$cvm) ### best lambda corresponds to min  cvm, mse.

    RSquare[[paste0("Factor_", i)]] <- r_squared(as.matrix(y), as.matrix(predict(cv.out, newx = train.data, s = cv.out$lambda.min)))

    CV_RMSE[[paste0("Factor_", i)]] <- sqrt(mse)
    lambda_out[[paste0("Factor_", i)]] <- cv.out$lambda.min

    Q22 <- lapply(unique(cv.out$foldid), function(id) {
      fit <- glmnet(x = train.data[cv.out$foldid != id,], ### this part is complicated, crossvalidiating the foldid, Fit excluding test set (foldid == id)
                    y = y[cv.out$foldid != id],
                    alpha = 0, standardize=FALSE)
      ## Test-set Y_hat using model fit at best lambda
      y_pred <- predict(fit, newx = train.data[cv.out$foldid == id,], s = cv.out$lambda.min)
      ## Test-set Y
      y_test <- y[cv.out$foldid == id]  ##  ## Test-set Y
      Q22 <- r_squared(y_test, y_pred) ## Cross-validated outer R2/ Q2
      return(Q22)})
    QSquare[[paste0("Factor_", i)]] <- mean(as.numeric(Q22)) ## Q22 lapply fun(id) calulates Q2 for ever fold, now we take the mean

    res <- do.call(cbind, list(RSquare, QSquare, CV_RMSE, lambda_out)) ## #bind all lists into results
    res <- as.data.frame(res)
    colnames(res) <- c('RSquared_Y','QSquared_Y', 'CV_RMSE', 'lambda_min')
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
    colnames(res) <- c('RSquared_Y', 'QSquared_Y', 'NComponents', 'CV_RMSEP', 'model')
  }
  return(res)
}

