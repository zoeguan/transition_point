source('transition_point_fns.R')
source('sim_fns.R')

library(caret)

# number of studies
ndat = 4
# size of each study
n = 40
# number of covariates
nvar = 10
# sigma_eps
sigma_eps = 1
# beta
set.seed(123)
beta = rnorm(nvar, 0, 1)
ind.small = sample(1:nvar, nvar/2)
beta[ind.small] = rnorm(length(ind.small), 0, 0.01)

load("edat_orig.RData")

# remove small studies 
edat_orig = edat_orig[which(sapply(edat_orig, function(x) nrow(x))>=50)]

# training set
set.seed(1)
edat_train <- init_data(edat_orig[1:6], ndat, nvar)
for (i in 1:length(edat_train)) {
  dataset = edat_train[[i]]
  edat_train[[i]] = dataset[1:n,]
}
# test set
set.seed(100)
edat_test <- init_data(edat_orig[7:13], ndat, nvar)
for (i in 1:length(edat_test)) {
  dataset = edat_test[[i]]
  edat_test[[i]] = dataset[1:n,]
}

### model hyperparameters 
# k-fold cross-validation for single study and merged dataset
set.seed(12)
study.num = sample(1:length(edat_train), 1)
train1 = edat_train[[study.num]]
eps1 = rnorm(nrow(train1), 0, sigma_eps)
train1$y = as.matrix(train1)[, 1:nvar] %*% beta + eps1
train2 = rbindlist(edat_train)
eps2 = rnorm(nrow(train2), 0, sigma_eps)
train2$y = as.matrix(train2)[, 1:nvar] %*% beta + eps2
# create folds
k=5
folds1 = createFolds(train1$y, k, list = TRUE, returnTrain = FALSE)
folds2 = createFolds(train2$y, k, list = TRUE, returnTrain = FALSE)
# LASSO, ridge
lambdas = 2^seq(-8, 8)
error.lasso = rep(0, length(lambdas))
error.ridge = rep(0, length(lambdas))
error.lasso.m = rep(0, length(lambdas))
error.ridge.m = rep(0, length(lambdas))
for (i in 1:length(lambdas)) {
  lambda = lambdas[i]
  for (j in 1:length(folds1)) {
    dataset = train1[unlist(folds1[-j]),]
    test = train1[folds1[[j]],]
    sd_y = sqrt(var(dataset$y)*(length(dataset$y)-1)/length(dataset$y))
    fit.lasso = glmnet(as.matrix(dataset)[, 1:nvar], dataset$y, alpha = 1, 
                       lambda = lambda,
                       intercept=F)
    fit.ridge = glmnet(as.matrix(dataset)[, 1:nvar], dataset$y, alpha = 0, 
                       lambda = lambda*sd_y/length(dataset$y),
                       intercept=F)
    error.lasso[i] = error.lasso[i] + mean((predict(fit.lasso, newx = as.matrix(test)[, 1:nvar]) - test$y)^2)
    error.ridge[i] = error.ridge[i] + mean((predict(fit.ridge, newx=as.matrix(test)[, 1:nvar]) - test$y)^2)
    
    dataset.m = train2[unlist(folds2[-j]),]
    test.m = train2[folds2[[j]],]
    sd_y = sqrt(var(dataset.m$y)*(length(dataset.m$y)-1)/length(dataset.m$y))
    fit.lasso = glmnet(as.matrix(dataset.m)[, 1:nvar], dataset.m$y, alpha = 1, 
                       lambda = lambda,
                       intercept=F)
    fit.ridge = glmnet(as.matrix(dataset.m)[, 1:nvar], dataset.m$y, alpha = 0, 
                       lambda = lambda*sd_y/length(dataset.m$y),
                       intercept=F)
    error.lasso.m[i] = error.lasso.m[i] + mean((predict(fit.lasso, newx = as.matrix(test.m)[, 1:nvar]) - test.m$y)^2)
    error.ridge.m[i] = error.ridge.m[i] + mean((predict(fit.ridge, newx=as.matrix(test.m)[, 1:nvar]) - test.m$y)^2)
  }
}
lambdak.lasso = lambdas[which(error.lasso==min(error.lasso))[1]]
lambdak.ridge = lambdas[which(error.ridge==min(error.ridge))[1]]
lambda.lasso = lambdas[which(error.lasso.m==min(error.lasso.m))[1]]
lambda.ridge = lambdas[which(error.ridge.m==min(error.ridge.m))[1]]
# random forest 
lm.formula2 = as.formula(paste("y~", paste(c(names(train1)[1:nvar]), collapse="+")))
mtries = seq(floor(nvar/3), nvar)
error.rf = rep(0, length(mtries))
error.rf.m = rep(0, length(mtries))
for (i in 1:length(mtries)) {
  mtry = mtries[i]
  for (j in 1:length(folds1)) {
    dataset = train1[unlist(folds1[-j]),]
    test = train1[folds1[[j]],]
    rf = randomForest(lm.formula2, data=dataset, mtry=mtry)
    error.rf[i] = error.rf[i] + mean((predict(rf, newdata=test) - test$y)^2)
    
    dataset.m = train2[unlist(folds2[-j]),]
    test.m = train2[folds2[[j]],]
    rf = randomForest(lm.formula2, data=dataset.m, mtry=mtry)
    error.rf.m[i] = error.rf.m[i] + mean((predict(rf, newdata=test.m) - test.m$y)^2)
  }

}
mtryk = mtries[which(error.rf==min(error.rf))[1]]
mtry = mtries[which(error.rf.m==min(error.rf.m))[1]]
# neural network 
sizes = c(5, 10, 15, 20)
decays = c(0:5)
error.nn = matrix(0, nrow=length(sizes), ncol=length(decays))
error.nn.m = matrix(0, nrow=length(sizes), ncol=length(decays))
for (i in 1:length(sizes)) {
  size = sizes[i]
  for (l in 1:length(decays)) {
    decay = decays[l]
    
    for (j in 1:length(folds1)) {
      dataset = train1[unlist(folds1[-j]),]
      test = train1[folds1[[j]],]
      fit.nn = nnet(lm.formula2, data = dataset, size = size, linout = T, trace = F, decay = decay)
      error.nn[i, l] = error.nn[i, l] + mean((as.vector(predict(fit.nn, newdata=test))-test$y)^2)
      
      dataset.m = train2[unlist(folds2[-j]),]
      test.m = train2[folds2[[j]],]
      fit.nn.m = nnet(lm.formula2, data = dataset.m, size = size, linout = T, trace = F, decay = decay)
      error.nn.m[i, l] = error.nn.m[i, l] + mean((as.vector(predict(fit.nn.m, newdata=test.m))-test.m$y)^2)
    }
  }
}
ind.best = which(error.nn == min(error.nn), arr.ind = TRUE)
sizek = sizes[ind.best[1]]
decayk = decays[ind.best[2]]
ind.best.m = which(error.nn.m == min(error.nn.m), arr.ind = TRUE)
size = sizes[ind.best.m[1]]
decay = decays[ind.best.m[2]]


save(ndat, n, nvar, sigma_eps, beta, 
     edat_train, edat_test,
     lambda.ridge, lambdak.ridge,
     lambda.lasso, lambdak.lasso,
     mtry, mtryk,
     size, sizek,
     decay, decayk,
     error.lasso, error.lasso.m,
     error.ridge, error.ridge.m,
     error.rf, error.rf.m,
     error.nn, error.nn.m,
     file="hyperparams_thm1.RData")
