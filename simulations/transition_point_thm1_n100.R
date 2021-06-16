source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1_n100.RData")

### transition point for least squares
# note: in the paper the scaling factor for the \overline{\sigma^2} is 1/P where P is the number of original predictors (ncol(edat_train_misspec[[1]])), but tau_range() uses 1/ncol(edat_train[[1]]), and here ncol(edat_train[[1]]) is the number of basis-expanded predictors, so sigma2_bar below should be multiplied by ncol(edat_train[[1]])/ncol(edat_train_misspec[[1]]) to be consistent with the paper
sigma2_bar = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=sigma_eps)
sigma2_star = (nvar/ncol(Z_train[[1]]))*sigma2_bar



### transition point for ridge
sigma2_bar.ridge = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, length(edat_train)), sigma_eps=sigma_eps)
sigma2_star.ridge = (nvar/ncol(Z_train[[1]]))*sigma2_bar.ridge


### run simulations
sigma.vals = seq(0, sigma2_star.ridge[1]*3, length=13)
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
  print(j)
  sigmas = rep(sqrt(sigma.vals[j]), ncol(Z_train[[1]]))

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


save(err, results, sigma.vals, sigma2_bar, sigma2_star,
     sigma2_bar.ridge, sigma2_star.ridge,
     file="transition_point_thm1_n100.RData")





