###################################
#### Simulation functions
###################################

library(data.table)
library(nlme)
library(glmnet)
library(randomForest)
library(nnet)
library(metafor)
library(doMC)
library(parallel)
#library(sandwich)
#library(lmtest)
library(matrixStats)
library(mvmeta)
library(MASS)

### function to sample datasets
# edat_orig: list of datasets
# ndat: number of datasets to sample
# nvar: number of predictors to sample
init_data <- function(edat_orig, ndat, nvar){
  
  edat <- edat_orig
  edat <- edat[sample(1:length(edat), ndat)] # Randomize dataset order
  
  idx <- sample(1:ncol(edat[[1]]), nvar)
  for(i in 1:ndat){
    edat[[i]] <- edat[[i]][,idx]
    edat[[i]] <- as.data.frame(edat[[i]])
    colnames(edat[[i]]) <- paste0("V", 1:nvar)
  }
  return(edat)
}

### function to calculate prediction error for various models for a given value of sigma_re
# edat_train: list of design matrices for training datasets
# edat_test: list of design matrices for test datasets 
# f_train: f(X_k), true mean of Y_k (X_k %*% beta in old simulations)
# f_test: f(X_0), true mean of Y_0 (X_0 %*% beta in old simulations)
# sigma_re: vector where entry i is the square root of the variance of the random effect for predictor i (if predictor i does not have a random effect, then set the corresponding entry to 0; can also provide a single number instead of the vector if all predictors have random effects with the same variance)
# sigma_eps: square root of variance of residual errors
# wk: vector of weights for CSL
# standardize: TRUE if the predictors should be standardized before running LASSO/ridge regression, FALSE otherwise
# lambda_lasso: LASSO regularization parameter for merged learner
# lambdak_lasso: vector of LASSO regularization parameters for single-study learners (can also provide a single number instead of a vector if this hyperparameter is the same for all studies)
# lambda_ridge: ridge regularization parameter for merged learner
# lambdak_ridge: vector of ridge regularization parameters for single-study learners (can also provide a single number instead of a vector if this hyperparameter is the same for all studies)
# size: neural network hidden layer size for merged learner
# sizek: vector of neural network hidden layer sizes for single-study learners (can also provide a single number instead of a vector if this hyperparameter is the same for all studies)
# decay: neural network regularization parameter for merged learner
# decayk: vector of neural network regularization parameters for single-study learners (can also provide a single number instead of a vector if the regularization parameter is the same for all studies)
# mtry: number of variables sampled at each split for merged random forest
# mtryk: vector of numbers of variables sampled at each split for single-study random forests (can also provide a single number instead of a vector if this hyperparameter is the same for all studies)
# me: whether or not to fit a mixed effects model
# ma: whether or not to run fixed effects and random effects meta-analyses
# ranfor: whether or not to fit a random forest
# nn: whether or not to fit a neural network
# ols: whether or not to fit OLS
# maxWeights: maximum number of weights for neural network
# correlated_re: 0 if random effects uncorrelated; -1 if negatively correlated and 1 if positively correlated
# opt_wk: whether to estimate optimal weights for least squares and ridge regression
sim_rep_np = function(edat_train, edat_test, f_train, f_test, Z_train, Z_test, sigma_re, sigma_eps, wk,
                   standardize=F, 
                   lambda_lasso=1, lambdak_lasso=1, 
                   lambda_ridge=1, lambdak_ridge=1, 
                   size=10, sizek=10, decay=2, decayk=2,
                   mtry=5, mtryk=5,
                   me=F, ma=F,
                   ranfor=T, nn=T,
                   ols=T, maxWeights=1000, correlated_re=0, opt_wk=F) {
  
  intercept = F
  intercept_re = F
  if (all(edat_train[[1]][,1]==1)) {
    intercept = T
    if (all(Z_train[[1]][,1]==1)) {
      intercept_re = T
    }
  }
  
  pred.err = data.frame(lme=NA, lme2=NA,
                        merged=NA, csl=NA,
                        lasso.m=NA, lasso.csl=NA, 
                        ridge.m=NA, ridge.csl=NA, 
                        rf.m=NA, rf.csl=NA, 
                        nn.m=NA, nn.csl=NA,
                        meta.re=NA, 
                        meta.multi.reml=NA)
  
  nvar = ncol(edat_train[[1]])
  
  ind.features = 1:nvar
  if (intercept) {
    ind.features = 2:nvar
  }
  
  # generate outcomes for training data
  # y_k = f(x_k) + z_k gamma_k + eps_k
  edat_sim = edat_train
  for (k in 1:length(edat_train)) {
    dataset = edat_train[[k]]
    f_k = as.matrix(f_train[[k]])
    Z_k = as.matrix(Z_train[[k]])
    gamma = rnorm(ncol(Z_k), 0, sigma_re)
    # make random effects correlated if needed
    if (correlated_re %in% 1) {
      sigma.mat = diag(sigma_re) + 0.25*min(sigma_re)
      diag(sigma.mat) = sigma_re
      gamma = mvrnorm(n=1, mu=rep(0, ncol(Z_k)), Sigma=sigma.mat)
    } else if (correlated_re %in% -1) {
      sigma.mat = - diag(sigma_re) - 0.25*min(sigma_re)
      diag(sigma.mat) = sigma_re
      gamma = mvrnorm(n=1, mu=rep(0, ncol(Z_k)), Sigma=sigma.mat)
    }
    eps = rnorm(nrow(dataset), 0, sigma_eps)
    dataset$y = f_k + 
      Z_k %*% gamma +
      eps 
    dataset$study = k
    edat_sim[[k]] = dataset
  }
  
  # generate outcomes for test data
  edat_sim_test = edat_test
  for (k in 1:length(edat_test)) {
    dataset = edat_test[[k]]
    f_0 = as.matrix(f_test[[k]])
    Z_0 = as.matrix(Z_test[[k]])
    gamma2 = rnorm(ncol(Z_k), 0, sigma_re)
    # make random effects correlated if needed
    if (correlated_re %in% 1) {
      sigma.mat = diag(sigma_re) + 0.25*min(sigma_re)
      diag(sigma.mat) = sigma_re
      gamma2 = mvrnorm(n=1, mu=rep(0, ncol(Z_k)), Sigma=sigma.mat)
    } else if (correlated_re %in% -1) {
      sigma.mat = - diag(sigma_re) - 0.25*min(sigma_re)
      diag(sigma.mat) = sigma_re
      gamma2 = mvrnorm(n=1, mu=rep(0, ncol(Z_k)), Sigma=sigma.mat)
    }
    eps2 = rnorm(nrow(dataset), 0, sigma_eps)
    dataset$y = f_0 + 
      Z_0 %*% gamma2 +
      eps2
    dataset$study = k
    edat_sim_test[[k]] = dataset
  }
  
  train = rbindlist(edat_sim)
  test = rbindlist(edat_sim_test)
  
  # formula for linear models
  lm.formula = as.formula(paste("y~", paste(c(names(train)[ind.features], "0"), collapse="+")))
  reStruct.formula = as.formula(paste(paste("~", paste(c(names(Z_train[[1]]), "0"), collapse="+")), "|study"))
  if (intercept) {
    lm.formula = as.formula(paste("y~", paste(c(names(train)[ind.features]), collapse="+")))
    reStruct.formula = as.formula(paste(paste("~", paste(c(names(Z_train[[1]])[-1], "0"), collapse="+")), "|study"))
    if (intercept_re) {
      reStruct.formula = as.formula(paste(paste("~", paste(c(names(Z_train[[1]])[-1]), collapse="+")),
                                           "|study"))
    }
  }
  
  # mixed effects model
  if (me) {
    # unequal variances
    fit.lme = tryCatch(do.call(lme, list(lm.formula, data=train,
                                         random=reStruct(reStruct.formula, pdClass="pdDiag"))), error = function(e) NA)
    # equal variances
    fit.lme2 = tryCatch(do.call(lme, list(lm.formula, data=train,
                            random=reStruct(reStruct.formula, pdClass="pdIdent"))),
                        error = function(e) NA)
    pred.lme = tryCatch(predict(fit.lme, newdata = test, level=0), error = function(e) rep(NA, nrow(test)))
    pred.err$lme[1] = mean((pred.lme-test$y)^2)
    pred.lme2 = tryCatch(predict(fit.lme2, newdata = test, level=0), error = function(e) rep(NA, nrow(test)))
    pred.err$lme2[1] = mean((pred.lme2-test$y)^2)
  }
  
  # merged OLS
  if (ols) {
    fit.merged = lm(lm.formula, data=train)
    pred.merged = predict(fit.merged, newdata = test)
    pred.err$merged[1] = mean((pred.merged-test$y)^2)
  }
  
  # CSL weights
  wk.ols = wk
  wk.ridge = wk
  
  # estimate optimal weights
  if (opt_wk) {
    f_train.pred = f_train
    for (i in 1:length(f_train)) {
      f_train.pred[[i]] = tryCatch(predict(fit.lme2, newdata = edat_train[[i]], level=0), error = function(e) rep(NA, nrow(edat_train[[i]])))
    }
    f_test.pred = f_test
    for (i in 1:length(f_test)) {
      f_test.pred[[i]] = tryCatch(predict(fit.lme2, newdata = edat_test[[i]], level=0), error = function(e) rep(NA, nrow(edat_test[[i]])))
    }
    wk.ols = optimal_weights_ridge(edat_train, edat_test, f_train.pred, f_test.pred, summary(fit.lme2)$sigma, sigma_re, 0)

    wk.ridge = optimal_weights_ridge(edat_train, edat_test, f_train.pred, f_test.pred, summary(fit.lme2)$sigma, sigma_re, lambdak_ridge)
    
  }
  
  # merged LASSO
  fit.lasso = glmnet(as.matrix(train)[, ind.features], train$y, alpha=1, 
                     lambda=lambda_lasso,
                     intercept=intercept,
                     standardize=standardize)
  pred.lasso = predict(fit.lasso, newx = as.matrix(test)[, ind.features])
  pred.err$lasso.m[1] = mean((pred.lasso - test$y)^2)
  
  # merged ridge
  # need to rescale lambda because of the formulation of the objective function used by glmnet
  sd_y = sqrt(var(train$y)*(length(train$y)-1)/length(train$y))
  fit.ridge = glmnet(as.matrix(train)[, ind.features], train$y, alpha = 0, 
                     lambda = lambda_ridge*sd_y/length(train$y),
                     intercept=intercept,
                     standardize=standardize)
  pred.ridge = predict(fit.ridge, newx=as.matrix(test)[, ind.features])
  pred.err$ridge.m[1] = mean((pred.ridge - test$y)^2)
  
  lm.formula2 = as.formula(paste("y~", paste(c(names(train)[ind.features]), collapse="+")))
  
  # merged random forest
  if (ranfor) {
    fit.rf = randomForest(lm.formula2, data=train, mtry=mtry)
    pred.rf = predict(fit.rf, newdata=test)
    pred.err$rf.m[1] = mean((pred.rf - test$y)^2)
  }
  
  # merged neural net
  if (nn) {
    fit.nn = nnet(lm.formula2, data = train, size = size, linout = T, trace = F, 
                  decay = decay, MaxNWts=maxWeights)
    pred.nn = as.vector(predict(fit.nn, newdata=test))
    pred.err$nn.m[1] = mean((pred.nn-test$y)^2)
  }
  
  
  ### cross-study learners
  ols.est = matrix(data=NA, nrow=length(edat_sim), ncol=nrow(test))
  lasso.est = matrix(data=NA, nrow=length(edat_sim), ncol=nrow(test))
  ridge.est = matrix(data=NA, nrow=length(edat_sim), ncol=nrow(test))
  rf.est = matrix(data=NA, nrow=length(edat_sim), ncol=nrow(test))
  nn.est = matrix(data=NA, nrow=length(edat_sim), ncol=nrow(test))
  
  beta.uni = matrix(data=NA, nrow=length(edat_sim), ncol=nvar)
  beta.uni.var = matrix(data=NA, nrow=length(edat_sim), ncol=nvar)
  beta.multi = matrix(data=NA, nrow=length(edat_sim), ncol=nvar)
  beta.multi.var = vector(mode="list", length=length(edat_sim))
  
  if (length(lambdak_lasso)==1) {
    lambdak_lasso = rep(lambdak_lasso, length(edat_sim))
  }
  if (length(lambdak_ridge)==1) {
    lambdak_ridge = rep(lambdak_ridge, length(edat_sim))
  }
  if (length(sizek)==1) {
    sizek = rep(sizek, length(edat_sim))
  }
  if (length(decayk)==1) {
    decayk = rep(decayk, length(edat_sim))
  }
  if (length(mtryk)==1) {
    mtryk = rep(mtryk, length(edat_sim))
  }
  
  for (i in 1:length(edat_sim)) {
    
    dataset = edat_sim[[i]]
    
    # OLS
    if (ols) {
      fit = lm(lm.formula, data=dataset)
      ols.est[i, ] = predict(fit, newdata = test)
    }
    
    # estimates for meta-analysis
    if (ma) {
      lm.uni = sapply(ind.features, function(x) summary(lm(as.formula(paste("y~", paste(c(names(train)[x], "0"), collapse="+"))), data=dataset))$coefficients[, 1:2])
      if (intercept) {
        lm.uni = cbind(summary(lm(as.formula("y~1"), data=dataset))$coefficients[, 1:2], lm.uni)
      }
      # estimates from univariate regressions
      beta.uni[i, ] = lm.uni[1, ]
      beta.uni.var[i, ] = (lm.uni[2, ])^2
      # estimates from multivariate regression
      beta.multi[i, ] = fit$coefficients
      beta.multi.var[[i]] = vcov(fit)
    }
    
    # lasso
    fit.lasso = glmnet(as.matrix(dataset)[, ind.features], dataset$y, alpha=1, 
                       lambda=lambdak_lasso[i],
                       intercept=intercept,
                       standardize=standardize)
    pred.lasso = predict(fit.lasso, newx = as.matrix(test)[, ind.features])
    lasso.est[i, ] = pred.lasso
    
    # ridge
    sd_y = sqrt(var(dataset$y)*(length(dataset$y)-1)/length(dataset$y))
    fit.ridge = glmnet(as.matrix(dataset)[, ind.features], dataset$y, alpha = 0,
                       lambda = lambdak_ridge[i]*sd_y/length(dataset$y),
                       intercept=intercept,
                       standardize=standardize)
    pred.ridge = predict(fit.ridge, newx=as.matrix(test)[, ind.features])
    ridge.est[i, ] = pred.ridge
    
    # random forest
    if (ranfor) {
      fit.rf = randomForest(lm.formula2, data=dataset, mtry=mtryk[i])
      pred.rf = predict(fit.rf, newdata=test)
      rf.est[i, ] = pred.rf
    }
    
    # nn
    if (nn) {
      fit.nn = nnet(lm.formula2, data=dataset, size=sizek[i], linout=T, trace=F, decay=decayk[i],
                    MaxNWts=maxWeights)
      pred.nn = as.vector(predict(fit.nn, newdata=test))
      nn.est[i, ] = pred.nn
    }
    
  }
  
  # combine predictions across studies 
  if (ols) {
    pred.err$csl[1] = mean((wk.ols %*% ols.est-test$y)^2)
  }
  pred.err$lasso.csl[1] = mean((colWeightedMeans(lasso.est, wk)-test$y)^2)
  pred.err$ridge.csl[1] = mean( (wk.ridge %*% ridge.est -test$y)^2)
  if (ranfor) {
    pred.err$rf.csl[1] = mean((colWeightedMeans(rf.est, wk)-test$y)^2)
  }
  if (nn) {
    pred.err$nn.csl[1] = mean((colWeightedMeans(nn.est, wk)-test$y)^2)
  }
  
  # meta-analyses
  if (ma) {
    
    # univariate
    beta.re = sapply(1:ncol(beta.uni), function(i) metafor::rma(beta.uni[,i], vi=beta.uni.var[,i], method="DL")$beta)
    pred.re = as.matrix(test)[, 1:nvar] %*% beta.re
    pred.err$meta.re[1] = mean((pred.re - test$y)^2)
    
    # multivariate
    meta.reml = tryCatch(mvmeta(beta.multi, S=beta.multi.var, method="reml", bscov="diag"), error = function(e) NA)
    pred.reml = tryCatch(as.matrix(test)[, 1:nvar] %*% t(meta.reml$coefficients), error = function(e) rep(NA, nrow(test)))
    pred.err$meta.multi.reml[1] = mean((pred.reml - test$y)^2)
    
  }
  
  return(pred.err)
  
}


### function to calculate prediction error for various models across multiple replicates, given a value of sigma_re 
# nreps: number of replicates
# n_cores: number of cores for parallel computation
sim_multi_np = function(edat_train, edat_test, f_train, f_test, Z_train, Z_test, sigma_re, sigma_eps, wk,
               nreps=100,
               standardize=F,
               lambda_lasso=1, lambdak_lasso=1, 
               lambda_ridge=1, lambdak_ridge=1, 
               size=10, sizek=10, decay=2, decayk=2,
               mtry=5, mtryk=5,
               me=F, ma=F,
               ranfor=T, nn=T,
               ols=T, maxWeights=1000,
               n_cores=1, correlated_re=0, opt_wk=F) {
  
  registerDoMC(cores=n_cores)
  
  results = foreach (j=1:nreps, .combine=rbind) %dopar% {
    sim_rep_np(edat_train=edat_train, edat_test=edat_test, f_train=f_train, f_test=f_test, Z_train=Z_train, Z_test=Z_test, 
               sigma_re=sigma_re, 
               sigma_eps=sigma_eps,
               wk=wk,
               standardize=standardize,
               lambda_lasso=lambda_lasso, lambdak_lasso=lambdak_lasso,
               lambda_ridge=lambda_ridge, lambdak_ridge=lambdak_ridge, 
               size=size, sizek=sizek,
               decay=decay, decayk=decayk,
               mtry=mtry, mtryk=mtryk,
               me=me, ma=ma, ranfor=ranfor, nn=nn, ols=ols, maxWeights=maxWeights,correlated_re=correlated_re, opt_wk=opt_wk)
    
    
  }
  return(results)
}








### function to bootstrap MSE ratio for two models
# nboot: number of bootstrap iterations
# results.df: dataframe formatted like the output of sim_multi()
# col1: column of results.df corresponding to first model
# col2: column of results.df corresponding to second model
boot_ci = function(nboot, results.df, col1, col2, seed=1) {
  set.seed(seed)
  results.df = results.df[which(results.df[,col1]>0 & results.df[, col2]>0),]
  avg.mse.merged = rep(NA, nboot)
  avg.mse.csl = rep(NA, nboot)
  for (i in 1:nboot) {
    boot.sample.ind = sample(1:nrow(results.df), nrow(results.df), replace=T)
    avg.mse.merged[i] = mean(results.df[boot.sample.ind, col1])
    avg.mse.csl[i] = mean(results.df[boot.sample.ind, col2])
  }
  return(list(avg.mse.merged, avg.mse.csl, avg.mse.csl/avg.mse.merged))
}

