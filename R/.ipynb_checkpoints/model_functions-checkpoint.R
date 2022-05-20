###rsquare function

r_squared <- function(y, yhat) {
    ybar <- mean(y)
    ## Total SS
    ss_tot <- sum((y - ybar)^2)
    ## Residual SS
    ss_res <- sum((y - yhat)^2)
    ## R^2 = 1 - ss_res/ ss_tot
   R2_Y <-  1 - (ss_res / ss_tot)
    return(R2_Y)
}


### RCV indexing function ### 
 ### find index with min diffrence  ### find index with min diffrence 

Ind.RCV.k <- function(RCV,k){
    Rcv.k.index <- which.minn(RCV, k)
    return(Rcv.k.index)
    }
                             

          
                             


### Ridge.k function 
Ridge.k <-function(Data.k, y, seed_val){ 
model_out <- list()
for (i in seq_along(Data.k)) {  
                            set.seed(seed_val)
                            train.data.k <- as.matrix(Data.k[[i]]) ### One of the training Data sets extrpolated 
                            lambdas <- 10^seq(3, -2, by = -.1) ### glmnet has some bult in lambda range. Here we have values ranging from λ=1010 to λ=10−2, essentially covering the full range of scenarios from the null model containing only the intercept, to the least squares fit.
                            cv.out = glmnet::cv.glmnet(x=train.data.k, y=as.matrix(y), lambda=lambdas, keep=TRUE, type.measure="mse", alpha = 0, standardize= FALSE)
                            model_out[[paste0("Factor_",i)]] <- cv.out### save each model iteritviely 
}
return(model_out)
    }

### PLS.k function 
PLS.k <- function(x, y, seed_val){ 
model_out <- list()

    for (i in seq_along(x)) {
    set.seed(seed_val)
    train.data <- x[[i]] ### One of the training Data sets extrpolated 
    model <- pls::plsr(as.matrix(y)~as.matrix(train.data), data=as.data.frame(train.data), validation = "CV", scale=FALSE) ## validation-cv, 10 fold cross val
    model_out[[paste0("Factor_",i)]] <- model ### save each model iteritviely 
   
}
return(model_out)
    }

### Ridge Function 
Ridge <-function(x, y, seed_val){ 
RSquare <- list()
CV_RMSEP <- list()
lambda_out <- list()
res <- list()
QSquare <- list()
for (i in seq_along(x)) {  
                            set.seed(seed_val)
                                train.data <- as.matrix(x[[i]]) ### One of the training Data sets extrpolated 
                                lambdas <- 10^seq(3, -2, by = -.1) ### glmnet has some bult in lambda range. Here we have values ranging from λ=1010 to λ=10−2, essentially covering the full range of scenarios from the null model containing only the intercept, to the least squares fit.
                                cv.out = glmnet::cv.glmnet(x=train.data, y=as.matrix(y), lambda=lambdas, keep=TRUE, type.measure="mse", alpha = 0, standardize= FALSE) # Fit ridge regression model on training data, default 10 fold crossvalidation
                                bestlam = cv.out$lambda.min  # Select lamda that minimizes training MSE
                                mse <- min(cv.out$cvm) ### best lambda corresponds to min  cvm, mse. 
    
                                RSquare[[paste0("Factor_", i)]] <- r_squared(as.matrix(y), as.matrix(predict(cv.out, newx = train.data, s = cv.out$lambda.min)))
                                
                                CV_RMSEP[[paste0("Factor_", i)]] <- sqrt(mse)
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
    
                                res <- do.call(cbind, list(RSquare, QSquare, CV_RMSEP, lambda_out)) ## #bind all lists into results 
                                res <- as.data.frame(res)
                                colnames(res) <- c('RSquared_Y','QSquared_Y', 'CV_RMSEP', 'lambda_min') 
}
        
return(res)
    }
### PLS Function 
PLS <- function(x, y, seed_val){ 
NComponents = list() ## define key outputs as lists first to save each for loop iteration
RSquare <- list() ## ^
QSquare <- list()
model_out <- list()
res <- list()
CV_RMSEP <- list()
    for (i in seq_along(x)) {
                    set.seed(seed_val)
                    train.data <- x[[i]] ### One of the training Data sets extrpolated 
                    model <- pls::plsr(as.matrix(y)~as.matrix(train.data), data=as.data.frame(train.data), ncomp=10, validation = "CV", segments = 10, scale=FALSE ) ## validation-cv, 10 fold cross val
                    CVRMSEP <- pls::RMSEP(model,'CV')$val ### cross-model-validated RMSEP calculated from the training data set- as not to leak info 
                    NComp <- which.min(CVRMSEP) ## index of minimun CVRMSEP to detemine optimal number of components 
                    RSquare[[paste0(names(x[i]))]] <- pls::R2(model, 'train')$val[as.numeric(NComp)] ## Take R2 for optimal no. comps  
                    QSquare[[paste0(names(x[i]))]]  <- pls::R2(model, 'CV')$val[as.numeric(NComp)] ## Q2 also known as cross-validated R2
                    CV_RMSEP[[paste0("Factor_",i)]]<- min(CVRMSEP)
                    NComponents[[paste0("Factor_",i)]]<- as.numeric(NComp)## save no.componenst 
                    res <- do.call(cbind, list(RSquare, QSquare, NComponents, CV_RMSEP)) ## #bind all lists into results
                    res <- as.data.frame(res)
                    colnames(res) <- c('RSquared_Y', 'QSquared_Y', 'NComponents', 'CV_RMSEP')
}
return(res)
    }

### Index min k values in array 

which.minn <- function(x,k){
  if (k==1)
    which.min(x[x>0])
  else
    {
      if (k>1){
        ii <- order(x,decreasing=FALSE)[1:min(k,length(x))]
        ii[!is.na(x[ii])]
      }
      else {
       stop("k must be >=1")
      }
    }
}
