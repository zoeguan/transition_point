library(gridExtra)
library(ggplot2)
library(reshape2)

load("transition_point_thm1_opt.RData")
source('transition_point_fns.R')
source('sim_fns.R')

ratio = data.frame(sigma2 = sigma2.vals,
                   csl=err$csl / err$merged,
                   ridge=err$ridge.csl / err$ridge.m)
ratio.lo = ratio
ratio.hi = ratio



### bootstrap CIs for MSE ratio comparing each CSL learner to the corresponding merged learner

for (i in 1:length(sigma2.vals)) {
  res.boot = boot_ci(200, results[[i]], "merged", "csl")
  ratio.lo$csl[i] = quantile(res.boot[[3]], 0.025)
  ratio.hi$csl[i] = quantile(res.boot[[3]], 0.975)
}

for (i in 1:length(sigma2.vals)) {
  res.boot = boot_ci(200, results[[i]], "ridge.m", "ridge.csl")
  ratio.lo$ridge[i] = quantile(res.boot[[3]], 0.025)
  ratio.hi$ridge[i] = quantile(res.boot[[3]], 0.975)
}



load("hyperparams_thm1.RData")

# scaling factor for sigma2_star, sigma2_star.ridge, etc to make the x-axis values consistent with the definition of \overline{\sigma^2} in paper
conversion.ratio = ncol(Z_train[[1]])/ncol(edat_train_misspec[[1]])

# scaling factor for sigma2_bar, sigma2_bar.ridge to make the x-axis values consistent with the definition of \overline{\sigma^2} in paper
conversion.ratio.2 = ncol(edat_train[[1]])/ncol(edat_train_misspec[[1]])

### transition point for least squares
sigma2_bar = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train)), sigma_eps=sigma_eps)*conversion.ratio.2

### transition point for least squares
sigma2_bar_ridge = tau_range(edat_train, edat_test, f_train, f_test, Z_train, lambda=lambda.ridge, lambdak=rep(lambdak.ridge, length(edat_train)), sigma_eps=sigma_eps)*conversion.ratio.2


sigma2.vals = sigma2.vals*conversion.ratio

### MSE ratio plots

bar.width = .003

p1 = qplot(sigma2.vals, log(ratio$csl), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_vline(xintercept=tpt_ls*conversion.ratio.2, color="#009900", linetype = "dotted", size=1.2) +
  geom_vline(xintercept=sigma2_bar, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$csl), ymax=log(ratio.hi$csl)), width=bar.width) +
  ggtitle("Least Squares") +
  theme(text = element_text(size=9))

p2 = qplot(sigma2.vals, log(ratio$ridge), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_vline(xintercept=tpt_r*conversion.ratio.2, color="#009900", linetype = "dotted", size=1.2) +
  geom_vline(xintercept=sigma2_bar_ridge, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$ridge), ymax=log(ratio.hi$ridge)), width=bar.width) +
  ggtitle("Ridge Regression") +
  theme(text = element_text(size=9))


png('plots_thm1_opt.png', width=600, height=200, res=100)
grid.arrange(p1, p2, ncol=2)
dev.off()
