########################################################
### Functions for calculating the transition point
########################################################

library(data.table)

### function to calculate matrix trace
tr = function(M) sum(diag(M))

### function to calculate transition bounds for ridge regression (Theorem 2)
# edat_train: tilde(X_k)
# edat_test: tilde(X_0)
# f_train: f(X_k)
# f_test: f(X_0)
# Z_train: random effects design matrix Z_k
# wk: vector of CSL weights 
# lambda: regularization parameter for merged learner (set to 0 for least squares)
# lambdak: vector of regularization parameters for single-study learners (set to 0 for least squares)
# sigma_eps: square root of variance of residual error terms
# cols_re_list: list of vectors of indices obtained by partitioning the predictors with random effects into subsets whose random effects have the same variance 
  # for example, if the predictors with random effects are 2, 3, and 4, and the random effects of 2 and 4 have the same variance while the random effect of 3 has a different variance than 2 and 4, then set cols_re_list=list(c(2, 4), 3)
tau_range = function(edat_train, edat_test, f_train, f_test, Z_train=edat_train, wk=rep(1/length(edat_train), length(edat_train)), lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=1, cols_re_list=list(1:ncol(Z_train[[1]]))) {
  
  # K = ndat
  ndat = length(edat_train)
  # P = nvar
  nvar = ncol(edat_train[[1]])
  
  # tilde(X_0), design matrix for test data
  X_0 = as.matrix(rbindlist(edat_test))
  
  f_0 = Reduce(c, f_test)
  f = Reduce(c, f_train)
  
  # tilde(X)
  X = as.matrix(rbindlist(edat_train))
  
  
  # I^-
  I_minus = diag(nvar)
  if (all(X[, 1]==1)) {
    I_minus[1, 1] = 0
  }
  
  # tilde(R) = tilde(X)^T tilde(X) + lambdak I^-
  R = t(X) %*% X + lambda * I_minus
  
  # tilde(R_k) = tilde(X_k)^T tilde(X_k) + lambdak I^-
  R_k.list = vector(mode="list", length=ndat)
  # weighted bias for study k = wk*(tilde(X_0) tilde(R_k)^-1 tilde(X_k)^T f(X_k) - f(X_0))
  bias_k.list = vector(mode="list", length=ndat)
  for (k in 1:ndat) {
    X_k = as.matrix(edat_train[[k]])
    f_k = as.matrix(f_train[[k]])
    R_k.list[[k]] = t(X_k) %*% X_k + lambdak[k] * I_minus
    bias_k.list[[k]] = wk[k]*(X_0 %*% solve(R_k.list[[k]]) %*% t(X_k) %*% f_k - f_0)
  }
  
  # bias terms
  # b_ens = sum_k w_k tilde(X_0) 
  b_ens = Reduce("+", bias_k.list)
  # b_merge = tilde(X_0) tilde(R)^{-1} tilde(X)^T f(X) - f(X_0)
  b_merge = X_0 %*% solve(R) %*% t(X) %*% f - f_0
  

  num1k = vector(mode="list", length=ndat)
  denom1k = vector(mode="list", length=ndat)
  denom2k = vector(mode="list", length=ndat)
  for (k in 1:ndat) {
    X_k = as.matrix(edat_train[[k]])
    Z_k = as.matrix(Z_train[[k]])
    R_k = R_k.list[[k]]
    # wk^2 * tr(X_0 R_k^{-1} X_k^T X_k X_0^T )
    num1k[[k]] = wk[k]^2 * tr(X_0 %*% solve(R_k) %*% t(X_k) %*% X_k %*% solve(R_k) %*% t(X_0))
    denom1k[[k]] = t(Z_k) %*% X_k %*% solve(R) %*% t(X_0) %*% X_0 %*% solve(R) %*% t(X_k) %*% Z_k
    denom2k[[k]] = wk[k]^2 * t(Z_k) %*% X_k %*% solve(R_k) %*% t(X_0) %*% X_0 %*% solve(R_k) %*% t(X_k) %*% Z_k
  }
  
  num1 = Reduce("+", num1k)
  denom1 = Reduce("+", denom1k)
  denom2 = Reduce("+", denom2k)
  
  num2 = tr(X_0 %*% solve(R) %*% t(X) %*% X %*% solve(R) %*% t(X_0))
  
  denom_mat = denom1 - denom2
  
  squared_bias_ens = t(b_ens) %*% b_ens
  squared_bias_m = t(b_merge) %*% b_merge
  
  diag_means = sapply(1:length(cols_re_list), function(i) mean(diag(denom_mat)[cols_re_list[[i]]]))
  
  return(c( (sigma_eps^2*(num1-num2) + squared_bias_ens - squared_bias_m)/(nvar*max(diag_means)),
            (sigma_eps^2*(num1-num2) + squared_bias_ens - squared_bias_m)/(nvar*min(diag_means))))

  
}


### function to calculate optimal CSL weights for ridge regression 
optimal_weights_ridge = function(edat_train, edat_test, f_train, f_test, sigma_eps, sigma_re, lambdak, scale_mat_k=NULL) {
  
  ndat = length(edat_train)
  nvar = ncol(edat_train[[1]])
  
  I = diag(nvar)
  if (all(edat_train[[1]][, 1]==1)) {
    I[1, 1] = 0
  }
  
  if (is.null(scale_mat_k)) {
    scale_mat_k = rep(list(diag(nvar)), ndat)
  }
  
  if (length(lambdak)==1) {
    lambdak = rep(lambdak, ndat)
  }
  
  if (length(sigma_re)==1) {
    sigma_re = rep(sigma_re, nvar)
  }
  
  Rk = vector(mode="list", length=ndat)
  Mk = vector(mode="list", length=ndat)
  vk = rep(NA, ndat)
  bk = vector(mode="list", length=ndat)
  
  test_mat = as.matrix(rbindlist(edat_test))
  
  for (k in 1:ndat) {
    Xk = as.matrix(edat_train[[k]])
    Rk[[k]] = t(Xk) %*% Xk
    Mk[[k]] = Rk[[k]] + lambdak[k]*I %*% solve(scale_mat_k[[k]]) %*% solve(scale_mat_k[[k]])
    
    vk[k] = tr(solve(Mk[[k]]) %*% Rk[[k]] %*% diag(sigma_re^2) %*% Rk[[k]] %*% solve(Mk[[k]]) %*%  t(test_mat) %*% test_mat) + 
      sigma_eps^2*tr(solve(Mk[[k]]) %*% Rk[[k]] %*% solve(Mk[[k]]) %*% t(test_mat) %*% test_mat)
    
    bk[[k]] = test_mat %*% solve(Mk[[k]]) %*% t(Xk) %*% as.matrix(f_train[[k]]) - Reduce(c, f_test)
      
  }
  
  C = matrix(NA, ncol=ndat,nrow=ndat)
  for (i in 1:ndat) {
    for (j in 1:ndat) {
      if(i==j) {
        C[i, j] = vk[i] + t(bk[[i]]) %*% bk[[i]]
      } else {
        C[i, j] =  t(bk[[i]] ) %*% bk[[j]] 
      }
    }
  }
  
  J = rep(1, ndat)
  
  return(c(solve(C) %*% J %*% solve(t(J) %*% solve(C) %*% J)))
}



##### old functions for linear generating model

### function to calculate least squares transition bounds under linear generating model
# edat_train: list of design matrices for training datasets
# edat_test: list of design matrices for test datasets 
# wk: vector of CSL weights 
# sigma_eps: square root of variance of residual error terms
# cols_re_list: list of vectors of indices obtained by partitioning the predictors with random effects into subsets whose random effects have the same variance 
# for example, if the predictors with random effects are 2, 3, and 4, and the random effects of 2 and 4 have the same variance while the random effect of 3 has a different variance than 2 and 4, then set cols_re_list=list(c(2, 4), 3)
tau_ls_range = function(edat_train, edat_test, wk, sigma_eps, cols_re_list=list(1:ncol(edat_train[[1]]))) {
  
  ndat = length(edat_train)
  nvar = ncol(edat_train[[1]])
  
  # X_0, design matrix for test data
  test_mat = as.matrix(rbindlist(edat_test))
  
  Rk = vector(mode="list", length=ndat)
  num1k = vector(mode="list", length=ndat)
  for (k in 1:ndat) {
    # X_k, design matrix for study k
    Xk = as.matrix(edat_train[[k]])
    # R_k = X_k^T X_k
    Rk[[k]] = t(Xk) %*% Xk
    # wk^2 * tr(R_k^{-1} X_0^T X_0)
    num1k[[k]] = wk[k]^2 * tr(test_mat %*% solve(Rk[[k]]) %*% t(test_mat))
  }
  
  Rk_sum = Reduce(`+`, Rk)
  num1 = Reduce(`+`, num1k)
  
  denom1_mat_k = vector(mode="list", length=ndat)
  for (i in 1:ndat) {
    # R_k R^{-1} X_0^T X_0 R^{-1} R_k
    denom1_mat_k[[i]] = Rk[[i]] %*% solve(Rk_sum) %*% t(test_mat) %*% test_mat %*% solve(Rk_sum) %*%
      Rk[[i]]
  }
  
  denom1_mat = Reduce(`+`, denom1_mat_k)
  num2 = tr(test_mat %*% solve(Rk_sum) %*% t(test_mat))
  denom2_mat = (t(test_mat) %*% test_mat)*sum(wk^2)
  denom_mat = denom1_mat - denom2_mat
  
  diag_means = sapply(1:length(cols_re_list), function(i) mean(diag(denom_mat)[cols_re_list[[i]]]))
  
  return(c((sigma_eps^2/nvar)*(num1-num2)/(max(diag_means)),
           (sigma_eps^2/nvar)*(num1-num2)/(min(diag_means))))
  
}

### function to calculate least squares transition point under linear generating model
# cols_re: vector of indices corresponding to predictors with random effects
tau_ls = function(edat_train, edat_test, wk, sigma_eps, cols_re=1:ncol(edat_train[[1]])) {
  return(tau_ls_range(edat_train, edat_test, wk, sigma_eps, cols_re_list=list(cols_re))[1])
}

### function to calculate ridge regression transition point under linear generating model
# lambda: regularization parameter for merged learner
# lambdak: vector of regularization parameters for single-study learners
# beta: vector of coefficients (use true coefficients if known, otherwise use estimates)
# scale_mat: scaling matrix for merged dataset applied prior to running ridge regression 
# scale_mat_k: list of scaling matrices for training datasets applied prior to running ridge regression
# if the intercept is included when running ridge regression, then the first column of the design matrix should be all 1's
tau_r_range = function(edat_train, edat_test, wk, sigma_eps, cols_re_list=list(1:ncol(edat_train[[1]])), 
                       lambda=1, lambdak=rep(1, length(edat_train)), beta, 
                       scale_mat=NULL, scale_mat_k=NULL) {
  
  ndat = length(edat_train)
  nvar = ncol(edat_train[[1]])
  
  # merged design matrix
  edat_train_all = as.matrix(rbindlist(edat_train))
  
  # I^-
  I = diag(nvar)
  if (all(edat_train_all[, 1]==1)) {
    I[1, 1] = 0
  }
  
  if (is.null(scale_mat)) {
    scale_mat = diag(nvar)
  }
  if (is.null(scale_mat_k)) {
    scale_mat_k = rep(list(diag(nvar)), ndat)
  }
  
  if (length(lambdak)==1) {
    lambdak = rep(lambdak, ndat)
  }
  
  test_mat = as.matrix(rbindlist(edat_test))
  
  Rk = vector(mode="list", length=ndat)
  # M_k = R_k + \lambdak * I^-
  Mk = vector(mode="list", length=ndat)
  Bk = vector(mode="list", length=ndat)
  denom2k = vector(mode="list", length=ndat)
  num1k = vector(mode="list", length=ndat)
  for (k in 1:ndat) {
    Xk = as.matrix(edat_train[[k]])
    Rk[[k]] = t(Xk) %*% Xk
    Mk[[k]] = Rk[[k]] + lambdak[k]*I %*% solve(scale_mat_k[[k]]) %*% solve(scale_mat_k[[k]])
    denom2k[[k]] = wk[k]^2 * Rk[[k]] %*% solve(Mk[[k]]) %*%  t(test_mat) %*% test_mat %*%
      solve(Mk[[k]]) %*% Rk[[k]]
    num1k[[k]] = wk[k]^2 * tr(test_mat %*% solve(Mk[[k]]) %*% Rk[[k]] %*% solve(Mk[[k]]) %*% t(test_mat))
    Bk[[k]] = wk[k]*solve(Mk[[k]]) %*% Rk[[k]]
  }
  
  Rk_sum = Reduce(`+`, Rk)
  num1 = Reduce(`+`, num1k)
  Bk_sum = Reduce(`+`, Bk)
  denom2 = Reduce(`+`, denom2k)
  
  M = Rk_sum + lambda*I %*% solve(scale_mat) %*% solve(scale_mat)
  
  denom1_k = vector(mode="list", length=ndat)
  for (i in 1:ndat) {
    denom1_k[[i]] = Rk[[i]] %*% solve(M) %*%  t(test_mat) %*% test_mat %*% solve(M) %*% Rk[[i]]
  }
  denom1 = Reduce(`+`, denom1_k)
  
  num2 = tr(test_mat %*% solve(M) %*% Rk_sum %*% solve(M) %*% t(test_mat))
  
  denom_mat = denom1 - denom2
  
  squared_bias_csl = beta %*% t(Bk_sum - diag(nvar)) %*% t(test_mat) %*% test_mat %*% (Bk_sum - diag(nvar)) %*% beta
  squared_bias_m = beta %*% t(solve(M) %*% Rk_sum - diag(nvar)) %*% t(test_mat) %*% test_mat %*% 
    (solve(M) %*% Rk_sum - diag(nvar)) %*% beta
  
  diag_means = sapply(1:length(cols_re_list), function(i) mean(diag(denom_mat)[cols_re_list[[i]]]))
  
  return(c( (sigma_eps^2*(num1-num2) + squared_bias_csl - squared_bias_m)/(nvar*max(diag_means)),
            (sigma_eps^2*(num1-num2) + squared_bias_csl - squared_bias_m)/(nvar*min(diag_means))))
}



### function to calculate transition point from Theorem 3
tau_r = function(edat_train, edat_test, wk, sigma_eps, cols_re=1:ncol(edat_train[[1]]), 
                 lambda=1, lambdak=rep(1, length(edat_train)), beta, 
                 scale_mat=NULL, scale_mat_k=NULL) {
  return(tau_r_range(edat_train, edat_test, wk, sigma_eps, cols_re_list=list(cols_re), 
                     lambda=lambda, lambdak=lambdak, beta, 
                     scale_mat=scale_mat, scale_mat_k=scale_mat_k)[1])
}

### function to calculate optimal CSL weights for LS under linear generating model
optimal_weights0 = function(edat_train, edat_test, sigma_eps, sigma_re) {
  
  nvar = ncol(edat_train[[1]])
  
  if (length(sigma_re)==1) {
    sigma_re = rep(sigma_re, nvar)
  }
  
  test_mat = as.matrix(rbindlist(edat_test))
  wk = rep(NA, length(edat_train))
  
  for (i in 1:length(edat_train)) {
    Xk = as.matrix(edat_train[[i]])
    wk[i] = 1/(tr(test_mat %*% diag(sigma_re^2) %*% t(test_mat)) + sigma_eps^2 * tr(test_mat %*% solve(t(Xk) %*% Xk) %*% t(test_mat)))
  }
  
  return(wk/sum(wk))
}

### function to calculate scaling matrices for training datasets if using glmnet with standardize=T for ridge regression
scale_k_glmnet = function(edat_train) {
  scale_mat_k = vector(mode="list", length=length(edat_train))
  for (i in 1:ndat) {
    scale_mat_k[[i]] = diag(apply(edat_train[[i]], 2, function(x) 1/sqrt(var(x)*(nrow(edat_train[[i]])-1)/nrow(edat_train[[i]]))))
    if (all(edat_train[[1]][,1]==1)) {
      scale_mat_k[[i]][1, 1] = 1
    }
  }
  return(scale_mat_k)
}

### function to calculate scaling matrix fo merged dataset if using glmnet with standardize=T for ridge regression
scale_glmnet = function(edat_train) {
  scale_mat = diag(apply(rbindlist(edat_train), 2, function(x) 1/sqrt(var(x)*(nrow(rbindlist(edat_train))-1)/nrow(rbindlist(edat_train)))))
  if (all(edat_train[[1]][,1]==1)) {
    scale_mat[1, 1] = 1
  }
  return(scale_mat)
}



### function to calculate optimal CSL weights for ridge regression under linear generating model
optimal_weights_ridge0 = function(edat_train, edat_test, sigma_eps, sigma_re, lambdak, beta, scale_mat_k=NULL) {
  
  ndat = length(edat_train)
  nvar = ncol(edat_train[[1]])
  
  I = diag(nvar)
  if (all(edat_train[[1]][, 1]==1)) {
    I[1, 1] = 0
  }
  
  if (is.null(scale_mat_k)) {
    scale_mat_k = rep(list(diag(nvar)), ndat)
  }
  
  if (length(lambdak)==1) {
    lambdak = rep(lambdak, ndat)
  }
  
  if (length(sigma_re)==1) {
    sigma_re = rep(sigma_re, nvar)
  }
  
  Rk = vector(mode="list", length=ndat)
  Mk = vector(mode="list", length=ndat)
  vk = rep(NA, ndat)
  bk = vector(mode="list", length=ndat)
  
  test_mat = as.matrix(rbindlist(edat_test))
  
  for (k in 1:ndat) {
    Xk = as.matrix(edat_train[[k]])
    Rk[[k]] = t(Xk) %*% Xk
    Mk[[k]] = Rk[[k]] + lambdak[k]*I %*% solve(scale_mat_k[[k]]) %*% solve(scale_mat_k[[k]])
    
    vk[k] = tr(solve(Mk[[k]]) %*% Rk[[k]] %*% diag(sigma_re^2) %*% Rk[[k]] %*% solve(Mk[[k]]) %*%  t(test_mat) %*% test_mat) + 
      sigma_eps^2*tr(solve(Mk[[k]]) %*% Rk[[k]] %*% solve(Mk[[k]]) %*% t(test_mat) %*% test_mat)
    
    bk[[k]] = lambdak[k]*test_mat %*% solve(Mk[[k]]) %*% I %*% solve(scale_mat_k[[k]]) %*% solve(scale_mat_k[[k]]) %*% beta
  }
  
  C = matrix(NA, ncol=ndat,nrow=ndat)
  for (i in 1:ndat) {
    for (j in 1:ndat) {
      if(i==j) {
        C[i, j] = vk[i] + t(bk[[i]]) %*% bk[[i]]
      } else {
        C[i, j] =  t(bk[[i]] ) %*% bk[[j]] 
      }
    }
  }
  
  J = rep(1, ndat)
  
  return(c(solve(C) %*% J %*% solve(t(J) %*% solve(C) %*% J)))
}

