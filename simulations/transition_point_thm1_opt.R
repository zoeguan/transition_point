source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")

### optimal weights transition point for least squares
# note: in the paper the scaling factor for the \overline{\sigma^2} is 1/P where P is the number of original predictors (ncol(edat_train_misspec[[1]])), but tau_range() uses 1/ncol(edat_train[[1]]), and here ncol(edat_train[[1]]) is the number of basis-expanded predictors, so sigma2_bar below should be multiplied by ncol(edat_train[[1]])/ncol(edat_train_misspec[[1]]) to be consistent with the paper
cols_re = 1:5
sigma2_bar = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=sigma_eps)[1]
sigma2_star = (nvar/length(cols_re))*sigma2_bar
sigma2_temp = rep(sigma2_star, nvar)
sigma2_temp[-cols_re] = 0
wk = optimal_weights_ridge(edat_train, edat_test, f_train, f_test, sigma_eps, sqrt(sigma2_temp), 0)
tpt_ls = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=sigma_eps, wk=wk)[1]

dec = 0.00001
tol = 0.00001
while (tpt_ls < sigma2_bar & abs(tpt_ls-sigma2_bar)>tol) {
  sigma2_bar = sigma2_bar-dec
  sigma2_star = (nvar/length(cols_re))*sigma2_bar
  sigma2_temp = rep(sigma2_star, nvar)
  sigma2_temp[-cols_re] = 0
  wk = optimal_weights_ridge(edat_train, edat_test, f_train, f_test, sigma_eps, sqrt(sigma2_temp), 0)
  tpt_ls = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=sigma_eps, wk=wk)[1]
}

sigma2_star # 0.1430314
tpt_ls # 0.05500394


### optimal weights transition point for ridge regression
sigma2_bar_r = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat), sigma_eps=sigma_eps)[1]
sigma2_star_r = (nvar/length(cols_re))*sigma2_bar_r
sigma2_temp_r = rep(sigma2_star_r, nvar)
sigma2_temp_r[-cols_re] = 0
wk = optimal_weights_ridge(edat_train, edat_test, f_train, f_test, sigma_eps, sqrt(sigma2_temp_r), lambdak.ridge)
tpt_r = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat), sigma_eps=sigma_eps, wk=wk)[1]

dec = 0.00001
tol = 0.00001
while (tpt_r < sigma2_bar_r & abs(tpt_r-sigma2_bar_r)>tol) {
  sigma2_bar_r = sigma2_bar_r-dec
  sigma2_star_r = (nvar/length(cols_re))*sigma2_bar_r
  sigma2_temp_r = rep(sigma2_star_r, nvar)
  sigma2_temp_r[-cols_re] = 0
  wk = optimal_weights_ridge(edat_train, edat_test, f_train, f_test, sigma_eps, sqrt(sigma2_temp_r), lambdak.ridge)
  tpt_r = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat), sigma_eps=sigma_eps, wk=wk)[1]
}

sigma2_star_r # 0.3320354
tpt_r # 0.1276975



### run simulations
sigma2.vals = seq(0, sigma2_star_r*3, length=13)
err = data.frame(lme=rep(NA, length(sigma2.vals)), lme2=NA,
                 merged=NA, csl=NA,
                 lasso.m=NA, lasso.csl=NA, 
                 ridge.m=NA, ridge.csl=NA, 
                 rf.m=NA, rf.csl=NA, 
                 nn.m=NA, nn.csl=NA,
                 meta.re=NA, 
                 meta.multi.reml=NA)
results = vector("list", length = length(sigma2.vals))

ind = 1:length(sigma2.vals)
nreps = 1000

set.seed(1)
for (j in ind) {
  print(j)
  sigmas = rep(0, nvar)
  sigmas[cols_re] = sqrt(sigma2.vals[j])
  results[[j]] = sim_multi_np(edat_train=edat_train, edat_test=edat_test, f_train=f_train, f_test=f_test, Z_train=Z_train, Z_test=Z_test, 
                              sigma_re=sigmas, 
                              sigma_eps=sigma_eps,
                              wk=rep(1/ndat, ndat),
                              standardize=F, 
                              nreps=nreps,
                              lambda_lasso=lambda.lasso, lambdak_lasso=lambdak.lasso,
                              lambda_ridge=lambda.ridge, lambdak_ridge=lambdak.ridge, 
                              mtry=mtry, mtryk=mtryk,
                              size=size, sizek=sizek,
                              decay=decay, decayk=decayk,
                              me=T, ma=F, ranfor=F, nn=F, n_cores=8, opt_wk=T)
  r = results[[j]]
  res = colMeans(r)
  err[j, ] = res
  print(err[j, ])
}

save(err, results, sigma2.vals, 
     sigma2_star_r, sigma2_star, tpt_ls, tpt_r, 
     file="transition_point_thm1_opt.RData")

