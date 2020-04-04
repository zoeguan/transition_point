source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")

### optimal weights transition point for least squares
set.seed(123)
cols_re = sample(1:nvar, 5)
sigma2_bar = tau_ls(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, cols_re=cols_re)
sigma2_star = (nvar/length(cols_re))*sigma2_bar
sigma2_temp = rep(sigma2_star, nvar)
sigma2_temp[-cols_re] = 0
wk = optimal_weights(edat_train, edat_test, sigma_eps, sqrt(sigma2_temp))
tpt_ls = tau_ls(edat_train, edat_test, wk, sigma_eps, cols_re=cols_re)

dec = 0.00001
tol = 0.00001
while (tpt_ls < sigma2_bar & abs(tpt_ls-sigma2_bar)>tol) {
  sigma2_bar = sigma2_bar-dec
  sigma2_star = (nvar/length(cols_re))*sigma2_bar
  sigma2_temp = rep(sigma2_star, nvar)
  sigma2_temp[-cols_re] = 0
  wk = optimal_weights(edat_train, edat_test, sigma_eps, sqrt(sigma2_temp))
  tpt_ls = tau_ls(edat_train, edat_test, wk, sigma_eps, cols_re=cols_re)
}

sigma2_star # 0.1002635
tpt_ls # 0.05012846


### optimal weights transition point for ridge regression
sigma2_bar_r = tau_r(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, 
                     cols_re=cols_re, beta=beta, 
                     lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat),
                     scale_mat = scale_glmnet(edat_train),
                     scale_mat_k = scale_k_glmnet(edat_train))
sigma2_star_r = (nvar/length(cols_re))*sigma2_bar_r
sigma2_temp_r = rep(sigma2_star_r, nvar)
sigma2_temp_r[-cols_re] = 0
wk = optimal_weights_ridge(edat_train, edat_test, sigma_eps, sqrt(sigma2_temp_r), rep(lambdak.ridge, ndat), beta, scale_mat_k=scale_k_glmnet(edat_train))
tpt_r = tau_r(edat_train, edat_test, wk, sigma_eps, 
                     cols_re=cols_re, beta=beta, 
                     lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat),
                     scale_mat = scale_glmnet(edat_train),
                     scale_mat_k = scale_k_glmnet(edat_train))

dec = 0.0001
tol = 0.0001
while (tpt_r < sigma2_bar_r & abs(tpt_r-sigma2_bar_r)>tol) {
  sigma2_bar_r = sigma2_bar_r-dec
  sigma2_star_r = (nvar/length(cols_re))*sigma2_bar_r
  sigma2_temp_r = rep(sigma2_star_r, nvar)
  sigma2_temp_r[-cols_re] = 0
  wk = optimal_weights_ridge(edat_train, edat_test, sigma_eps, sqrt(sigma2_temp_r), rep(lambdak.ridge, ndat), beta, scale_mat_k=scale_k_glmnet(edat_train))
  tpt_r = tau_r(edat_train, edat_test, wk, sigma_eps, 
                       cols_re=cols_re, beta=beta, 
                       lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat),
                       scale_mat = scale_glmnet(edat_train),
                       scale_mat_k = scale_k_glmnet(edat_train))
}

sigma2_star_r # 0.09543735
tpt_r # 0.04770217



### run simulations
sigma2.vals = seq(0, sigma2_star*3, length=13)
err = data.frame(lme2=rep(NA, length(sigma2.vals)),
                 merged=NA, csl=NA,
                 ridge.m=NA, ridge.csl=NA)
results = vector("list", length = length(sigma2.vals))

ind = 1:length(sigma2.vals)
nreps = 1000

set.seed(1)
for (j in ind) {
  print(j)
  sigmas = rep(0, nvar)
  sigmas[cols_re] = sqrt(sigma2.vals[j])
  results[[j]] = sim_multi_weights(edat_train=edat_train, edat_test=edat_test, beta=beta, 
                                   wk=rep(1/ndat, ndat),
                                   standardize=T, 
                                   sigma_re=sigmas, 
                                   sigma_eps=sigma_eps,
                                   nreps=nreps,
                                   lambda_ridge=lambda.ridge, lambdak_ridge=lambdak.ridge, 
                                   n_cores=1, opt_wk=T)
  r = results[[j]]
  res = colMeans(r)
  err[j, ] = res
}

save(err, results, cols_re, sigma2.vals, 
     sigma2_star_r, sigma2_star, tpt_ls, tpt_r, 
     file="transition_point_thm1_opt.RData")

