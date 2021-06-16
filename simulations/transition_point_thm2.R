source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")

### transition bounds
# note: in the paper the scaling factor for the \overline{\sigma^2} is 1/P where P is the number of original predictors (ncol(edat_train_misspec[[1]])), but tau_range() uses 1/ncol(edat_train[[1]]), and here ncol(edat_train[[1]]) is the number of basis-expanded predictors, so sigma2_bar below should be multiplied by ncol(edat_train[[1]])/ncol(edat_train_misspec[[1]]) to be consistent with the paper
cols_re = 1:5

sigma2_bar = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train_misspec)), sigma_eps=sigma_eps, cols_re_list=as.list(cols_re))

sigma2_star_lo = (nvar/ncol(Z_train[[1]]))*sigma2_bar[1]
sigma2_star_hi = (nvar/ncol(Z_train[[1]]))*sigma2_bar[2]


sigma2_bar_ridge = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, length(edat_train_misspec)), sigma_eps=sigma_eps, cols_re_list=as.list(cols_re))
sigma2_bar_ridge
sigma2_star_ridge_lo = (nvar/ncol(Z_train[[1]]))*sigma2_bar_ridge[1]
sigma2_star_ridge_hi = (nvar/ncol(Z_train[[1]]))*sigma2_bar_ridge[2]


### run simulations
sigma.vals = seq(0, sigma2_star_ridge_hi*2, length=13) 
err = data.frame(lme=rep(NA, length(sigma.vals)), lme2=NA,
                 merged=NA, csl=NA,
                 lasso.m=NA, lasso.csl=NA, 
                 ridge.m=NA, ridge.csl=NA, 
                 rf.m=NA, rf.csl=NA, 
                 nn.m=NA, nn.csl=NA,
                 meta.re=NA, 
                 meta.multi.reml=NA)
results = vector("list", length = length(sigma.vals))



set.seed(100)
ind = 1:length(sigma.vals)
nreps = 1000
for (j in ind) {
  sigmas = sqrt(sigma.vals[j]*c(0.5, 0.5, 1, 1, 2))
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
                              me=F, ma=F, ranfor=T, nn=T, n_cores=8)
  err[j, ] = colMeans(results[[j]])
  print(err[j, ])
}

save(err, results, sigma.vals, sigma2_bar, sigma2_bar_ridge, cols_re,
     file="transition_point_thm2.RData")





