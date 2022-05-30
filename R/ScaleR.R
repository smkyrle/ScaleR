#' @title  ScaleR: An R package for investigating scaling factors effects on models
#'
#' @description This package employs GLM regression or PLS over multiple scaling factors. Default is GLM
#'
#' @param x is a data matrix of predictors, features in columns and samples in rows.
#'
#' @parma y is the outcome variable.
#'
#' @param model is the model type you want to employ (model='glm'), default is GLM.
#'
#' @param gamma is the Scaling factor parameter, e.g list(0, 0.5, 1), default is sequence between 0 and 1 with an interval of 0.1
#'
#' @param lambda is the lambda parameter for glmnet
#'
#' @param alpha is the alpha parameter for glmnet. Input 0 for ridge regression, 1 for lasso regression, and a sequence of numbers between 0 and 1 for elastic nets.
#'
#' @param ncomp is the number of components used in pls
#'
#' @param nfolds is the number of folds used in k-fold cross validation
#'
#' @param rounds is the number of rounds of gamma sampling. After each round the grid space either side of the model with the lowest MSE is searched in the next round. Default rounds = 1.
#'
#' @return dataframe of parameters and performance metrics (MSE, R2, Q2)
#'
#' @export

library(glmnet)

ScaleR.scale <- function(train, test, gamma){
    # Scale train and test data with input gamma value

    train.mean <- sapply(train, mean)
    train.sd <- sapply(train, sd)
    indicies <- which(train.sd != 0)

    train.sd <- train.sd[indicies]
    train.mean <- train.mean[indicies]
    train <- train[,indicies]
    test <- test[,indicies]
    train.scaled <- mapply(function(x,y,z) (x-y)/(z^gamma), train, train.mean, train.sd)
    test.scaled <- mapply(function(x,y,z) (x-y)/(z^gamma), test, train.mean, train.sd)
    return(list(train=train.scaled,test=test.scaled))
}

kfold <- function(x, k){
    # kfold cross val split function
    # Function for creating random kfold splits
    # Input dataset and number of folds
    # Returns fold indexes

    # Get number of samples and turn into list of indexes
    num_samples <- dim(x)[1]
    indexes = c(1:num_samples)
    remainder <- num_samples %% k # Left over number of samples from even split

    # List of number of samples in each split
    split_size <- c()
    for(i in 1:k) {
        if(remainder != 0){
            split_size[i] <- as.integer(num_samples/k) + 1
            remainder <- remainder - 1
        }
        else{
            split_size[i] <- as.integer(num_samples/k)
        }
    }

    # Randomly sample 5 splits - size of each split determined by split size list
    splits <- split(indexes,sample(rep(1:k,times=split_size)))

    # Generate k train-test pairs from splits
    folds <- vector('list', k)
    for(f in 1:k){
        test_split <- c(unlist(splits[f], recursive=FALSE), use.names=FALSE)
        train <- c(unlist(splits[-f], recursive=FALSE), use.names=FALSE)
        folds[[f]] <- list(train=train, test=test_split)
    }

    # Return kfold splits as list of train and test indexes
    return(folds)
}

fit.linear.model <- function(x, y, params){
    # Wrapper for glmnet
    lambda <- params$lambda
    alpha <- params$grid$alpha
    family <- params$family
    model <- glmnet::glmnet(x=as.matrix(x), y=as.matrix(y),lambda=lambda, alpha=alpha, standardize=FALSE, family=family)
    return(model)
}

fit.pls.model <- function(x, y, params){
    # Wrapper for plsr
    ncomp <- params$grid$ncomp
    model <- pls::plsr(as.matrix(y)~as.matrix(x), data=as.data.frame(x), ncomp=ncomp, validation = "none", scale=FALSE)
    return(model)
}

predict.model <- function(model, x, params){
    # Predict y for given x
    if(exists('lambda', where=params)){
        preds <- predict(model, newx = x)
    }
    else{
        preds <- predict(model, x, ncomp=params$grid$ncomp)
    }
    return(preds)
}

model.mse <- function(preds, y){
    # Calculate mean squared error
    res <- assess.glmnet(preds, newy=y)
    mse <- res$mse
    return(mse)
}

residuals.squared <- function(y, yhat, ybar) {
    # Calculate R2 or Q2 for given input
    ## Total SS
    ss_tot <- sum((y - ybar)^2)
    ## Residual SS
    ss_res <- sum((y - yhat)^2)
    ## R^2 = 1 - ss_res/ ss_tot
    R2_Y <-  1 - (ss_res / ss_tot)
    return(R2_Y)
}

fold.means <- function(y, folds){
    # calculate mean label value of training sets for R2 and Q2 calculations
    means<-list()
    for(fold in 1:length(folds)){
        fold <- folds[[fold]]
        mean <- mean(y[fold$train])
        means <- c(mean, means)
    }
    return(means)
}

prepare.grid <- function(gamma=NaN, alpha=NaN, lambda=NaN, ncomp=NaN){
    grid <- expand.grid(gamma=gamma, alpha=alpha, lambda=lambda, ncomp=ncomp)
    grid <- data.frame(t(na.omit(t(grid))))
    return(grid)
}

prepare.dataframe <- function(grid, nfolds){
    # Create empty dataframe for results

    ncol <- dim(grid)[2] + 11
    nrow <- nrow(grid)*nfolds
    df <- data.frame(matrix(ncol=ncol, nrow=nrow))
    colnames(df) <- c('index', 'fold', colnames(grid), 'mse', 'mse.mean', 'mse.sd','r2', 'r2.mean', 'r2.sd', 'q2', 'q2.mean', 'q2.sd')

    # Populate with fold and index data
    df$index <- rep(1:nrow(grid), each=nfolds)
    df$fold <- rep(1:nfolds, nrow(grid))

    for(header in colnames(grid)){
        column <- grid[[header]]
        df[header] <- rep(column, each=nfolds)
    }

    return(df)
}

process.results <-function(results){
    # Calculate mean and sd of mse, r2, and q2 across folds
    start <- 1
    nfolds <- max(results$fold)
    for(row in 1:(nrow(results)/nfolds)){
        end <- start - 1 + nfolds
        results$mse.mean[start:end] <- rep(mean(results$mse[start:end]),nfolds)
        results$mse.sd[start:end] <- rep(sd(results$mse[start:end]),nfolds)

        results$r2.mean[start:end] <- rep(mean(results$r2[start:end]),nfolds)
        results$r2.sd[start:end] <- rep(sd(results$r2[start:end]),nfolds)

        results$q2.mean[start:end] <- rep(mean(results$q2[start:end]),nfolds)
        results$q2.sd[start:end] <- rep(sd(results$q2[start:end]),nfolds)

        start <- end+1
    }
    return(results)
}

grid.search <- function(x, y, model, grid, params, folds){

    nfolds <- length(folds)

    # Prepare results df
    results <- prepare.dataframe(grid, nfolds)

    # Reset grid for grid search

    if(exists('lambda', where=params)){grid <- subset(grid, select = -lambda)}
    grid <- grid[!duplicated(grid), ]
    means <- fold.means(y, folds)
    perc_new <- NaN
    perc_old <- 0
    # Loop though grid rows

    for(i in 1:nrow(grid)){
        perc_new <- i/nrow(grid)*100
        if(perc_new-5 > perc_old){
            print(perc_new)
            perc_old <- perc_new
        }
        # Assign grid row to params
        params$grid <- grid[i,]

        # Loop over kfolds
        for(j in 1:nfolds){

            # select train and test folds
            # scale x data
            fold <- folds[[j]]
            x_tmp <- list(train=x[fold$train,], test=x[fold$test,])
            x_tmp <- ScaleR.scale(x_tmp$train, x_tmp$test, params$grid$gamma)
            y_tmp <- list(train=y[fold$train],test=y[fold$test])

            # fit chosen model
            fitted.model <- model(x_tmp$train, y_tmp$train, params=params)

            # test model
            test.preds <- predict.model(fitted.model, x_tmp$test, params)
            train.preds <- predict.model(fitted.model, x_tmp$train, params)
            mse <- model.mse(test.preds, y_tmp$test)

            # add mse to results dataframe
            if('lambda' %in% names(params)){
                for(l in 1:length(params$lambda)){

                    q2 <- residuals.squared(y_tmp$test, test.preds[,l], means[[j]])
                    r2 <- residuals.squared(y_tmp$train, train.preds[,l],  means[[j]])
                    results['mse'][results['lambda'] == params$lambda[l] & results['alpha'] == params$grid$alpha & results['gamma'] == params$grid$gamma & results['fold'] == j] = mse[l]
                    results['q2'][results['lambda'] == params$lambda[l] & results['alpha'] == params$grid$alpha & results['gamma'] == params$grid$gamma & results['fold'] == j] = q2
                    results['r2'][results['lambda'] == params$lambda[l] & results['alpha'] == params$grid$alpha & results['gamma'] == params$grid$gamma & results['fold'] == j] = r2
                }
            }
            else{
                q2 <- residuals.squared(y_tmp$test, test.preds, means[[j]])
                r2 <- residuals.squared(y_tmp$train, train.preds,  means[[j]])
                results['mse'][results['gamma'] == params$grid$gamma & results['fold'] == j & results['ncomp'] == params$grid$ncomp] = mse
                results['q2'][results['gamma'] == params$grid$gamma & results['fold'] == j & results['ncomp'] == params$grid$ncomp] = q2
                results['r2'][results['gamma'] == params$grid$gamma & results['fold'] == j & results['ncomp'] == params$grid$ncomp] = r2
            }

        }

    }

    # Calculate cv mean and sd for mse, r2, q2
    results <- process.results(results)

    return(results)
}

ScaleR <- function(x, y, model='glm', gamma=seq(0, 1, by=0.1), lambda=NULL, alpha=NULL, ncomp=NULL, nfolds=2, rounds=1){
    # params list for non-specific model input
    params <-list()

    # Prepare model
    if(model == 'glm'){
        model <- fit.linear.model
        ncomp <- NaN
        if(is.null(lambda)){
            lambda <- 10^seq(3, -2, by = -.1)
        }
        if(is.null(alpha)){
            alpha <- 0

        }
        params$lambda <- lambda
        print(alpha)
        if(length(alpha)>1){
            print("Running elastic nets")
        }
        else if(alpha==0){
            print("Running ridge regression")
        }
        else if(alpha==1){
            print("Running lasso regression")
        }

        if(length(unique(y)) == 2){
            cat('Two unique values identified. \nSetting family to binomial.\n')
            params$family = 'binomial'
        } else if(length(unique(y)) > 2){
            cat('More than two unique values identified. \nSetting family to gaussian.\n')
            params$family = 'gaussian'
        }
    }
    else if(model == 'pls'){
        model <- fit.pls.model
        if(is.null(ncomp)){
            ncomp <- c(1:20)
        }
        lambda <- NaN
        alpha <- NaN
    }

    # setup folds for cross val
    folds <- kfold(x, k=nfolds)
    results <- data.frame()
    interval <- gamma[2] - gamma[1]
    if(length(gamma) %% 2 == 0){
        len <- length(gamma)
    }
    else{
        len <- length(gamma) - 1
    }

    for(round in 1:rounds){
        grid <- prepare.grid(gamma=gamma, alpha=alpha, lambda=lambda, ncomp=ncomp)
        round_results <- grid.search(x, y, model=model, grid=grid, params=params, folds=folds)
        if(nrow(results)  != 0){
            round_results$index <- round_results$index + max(results$index)
        }
        results <- rbind(results, round_results)
        min.mse.mean <- round_results[which.min(round_results$mse.mean),]
        max.gamma <- min.mse.mean$gamma + interval/2
        min.gamma <- min.mse.mean$gamma - interval/2
        interval <- (max.gamma - min.gamma)/(len)
        #gamma = seq(min.gamma, max.gamma, by=interval)
    }
    return(results)
}
