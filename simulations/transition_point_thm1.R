source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")

### transition point
set.seed(123)
cols_re = sample(1:nvar, 5)
sigma2_bar = tau_ls(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, cols_re=cols_re)
sigma2_star = (nvar/length(cols_re))*sigma2_bar
sigma2_star
# 0.1045435
sigma2_bar
# 0.05227173

sigma2_bar_ridge = tau_r(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, 
                         cols_re=cols_re, beta=beta, 
                         lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat),
                         scale_mat = scale_glmnet(edat_train),
                         scale_mat_k = scale_k_glmnet(edat_train))
(nvar/length(cols_re))*sigma2_bar_ridge
# 0.1020373
sigma2_bar_ridge 
# 0.05101867


### run simulations
sigma.vals = seq(0, sigma2_star*3, length=13)
err = data.frame(lme=rep(NA, length(sigma.vals)), lme2=NA,
                 merged=NA, csl=NA,
                 lasso.m=NA, lasso.csl=NA, 
                 ridge.m=NA, ridge.csl=NA, 
                 rf.m=NA, rf.csl=NA, 
                 nn.m=NA, nn.csl=NA,
                 meta.re=NA, 
                 meta.multi.reml=NA)
results = vector("list", length = length(sigma.vals))

set.seed(12)
ind = 1:length(sigma.vals)
nreps = 1000
for (j in ind) {
  sigmas = rep(0, nvar)
  sigmas[cols_re] = sqrt(sigma.vals[j])
  results[[j]] = sim_multi(edat_train=edat_train, edat_test=edat_test, beta=beta, 
                           wk=rep(1/ndat, ndat),
                           standardize=T, 
                           sigma_re=sigmas, 
                           sigma_eps=sigma_eps,
                           nreps=nreps,
                           lambda_lasso=lambda.lasso, lambdak_lasso=lambdak.lasso,
                           lambda_ridge=lambda.ridge, lambdak_ridge=lambdak.ridge, 
                           mtry=mtry, mtryk=mtryk,
                           size=size, sizek=sizek,
                           decay=decay, decayk=decayk,
                           me=T, ma=T)
  err[j, ] = colMeans(results[[j]])
  print(err[j, ])
}


save(err, results, sigma.vals, sigma2_star, cols_re,
     sigma2_bar, sigma2_bar_ridge,
     file="transition_point_thm1.RData")





