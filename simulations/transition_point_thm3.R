source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm3.RData")

# this scenario is time-consuming so we include the option of running simulations for one sigma value at a time (instead of all 13 values) by using the command line to specify the index of the value of interest
ind = 1:13
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  ind = as.numeric(args[1])
}

### transition point
set.seed(100)
cols_re = sample(1:nvar, 10)
sigma2_bar = tau_r(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, 
                         cols_re=cols_re, beta=beta, 
                         lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat),
                         scale_mat = scale_glmnet(edat_train),
                         scale_mat_k = scale_k_glmnet(edat_train))
sigma2_bar
# 0.06109191
sigma2_star = (nvar/length(cols_re))*sigma2_bar
sigma2_star
# 0.6109191

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

set.seed(10)
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
                           me=F, ma=F,
                           ols=F, maxWeights=3100)
  r = results[[j]]
  res = colMeans(r)
  err[j, ] = res
  print(res)
}

save(err, results, 
     sigma.vals, sigma2_star, cols_re,
     sigma2_bar,
     file=paste0("thm3_thm4_RData/transition_point_thm3_", ind[1], "_", ind[length(ind)], ".RData"))






