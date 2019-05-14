source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")

### transition point
set.seed(123)
cols_re = sample(1:nvar, 5)
sigma2_star = sqrt(rnorm(nvar, 0.05, 0.5)^2)
sigma2_star[-cols_re] = 0
sigma2_bar = tau_ls_range(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, cols_re_list=as.list(cols_re))
sigma2_bar
sigma2_star_lo = nvar*sigma2_bar[1]*sigma2_star / sum(sigma2_star)
sigma2_star_hi = nvar*sigma2_bar[2]*sigma2_star / sum(sigma2_star)


sigma2_bar_ridge = tau_r_range(edat_train, edat_test, rep(1/ndat, ndat), sigma_eps, 
                               cols_re_list=as.list(cols_re), beta=beta, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, ndat), scale_mat = scale_glmnet(edat_train),
                               scale_mat_k = scale_k_glmnet(edat_train))
sigma2_bar_ridge 
# 0.04165740 0.07054708

### run simulations
scale.vals.2 = seq(nvar*sigma2_bar[1]/sum(sigma2_star), nvar*sigma2_bar[2]/sum(sigma2_star), length=4)
scale.vals.1 = rev(seq(nvar*sigma2_bar[1]/sum(sigma2_star), 0, by=scale.vals.2[1]-scale.vals.2[2]))
scale.vals.3 = seq(nvar*sigma2_bar[2]/sum(sigma2_star), 2*nvar*sigma2_bar[2]/sum(sigma2_star), 
                   by=scale.vals.2[2]-scale.vals.2[1])
scale.vals = c(scale.vals.1[1:(length(scale.vals.1)-1)], scale.vals.2, scale.vals.3[2:length(scale.vals.3)])[1:13]
err = data.frame(lme=rep(NA, length(scale.vals)), lme2=NA,
                 merged=NA, csl=NA,
                 lasso.m=NA, lasso.csl=NA, 
                 ridge.m=NA, ridge.csl=NA, 
                 rf.m=NA, rf.csl=NA, 
                 nn.m=NA, nn.csl=NA,
                 meta.re=NA, 
                 meta.multi.reml=NA)
results = vector("list", length = length(scale.vals))

set.seed(100)
ind = 1:length(scale.vals)
nreps = 1000
for (j in ind) {
  sigmas = sqrt(scale.vals[j]*sigma2_star)
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

save(err, results, scale.vals, sigma2_star, sigma2_bar, sigma2_bar_ridge,
     scale.vals.2, cols_re,
     sigma2_star_lo, sigma2_star_hi,
     file="transition_point_thm2.RData")





