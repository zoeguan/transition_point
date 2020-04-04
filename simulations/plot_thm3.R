library(gridExtra)
library(ggplot2)
library(reshape2)

source('transition_point_fns.R')
source('sim_fns.R')
load("transition_point_thm3.RData")
load("hyperparams_thm3.RData")

ratio = data.frame(sigma2 = sigma.vals,
                   csl=err$csl / err$merged,
                   lasso=err$lasso.csl / err$lasso.m,
                   ridge=err$ridge.csl / err$ridge.m,
                   rf=err$rf.csl / err$rf.m,
                   nn=err$nn.csl / err$nn.m)
ratio.lo = ratio
ratio.hi = ratio


### bootstrap CIs for MSE ratio comparing each CSL learner to the corresponding merged learner

for (i in 1:length(sigma.vals)) {
  res.boot = boot_ci(200, (results[[i]]), "lasso.m", "lasso.csl")
  ratio.lo$lasso[i] = quantile(res.boot[[3]], 0.025, na.rm=T)
  ratio.hi$lasso[i] = quantile(res.boot[[3]], 0.975, na.rm=T)
}

for (i in 1:length(sigma.vals)) {
  res.boot = boot_ci(200, (results[[i]]), "ridge.m", "ridge.csl")
  ratio.lo$ridge[i] = quantile(res.boot[[3]], 0.025, na.rm=T)
  ratio.hi$ridge[i] = quantile(res.boot[[3]], 0.975, na.rm=T)
}

for (i in 1:length(sigma.vals)) {
  res.boot = boot_ci(200, (results[[i]]), "nn.m", "nn.csl")
  ratio.lo$nn[i] = quantile(res.boot[[3]], 0.025, na.rm=T)
  ratio.hi$nn[i] = quantile(res.boot[[3]], 0.975, na.rm=T)
}

for (i in 1:length(sigma.vals)) {
  res.boot = boot_ci(200, (results[[i]]), "rf.m", "rf.csl")
  ratio.lo$rf[i] = quantile(res.boot[[3]], 0.025, na.rm=T)
  ratio.hi$rf[i] = quantile(res.boot[[3]], 0.975, na.rm=T)
}


### MSE ratio plots

sigma.vals = sigma.vals*length(cols_re)/nvar
sigma2_star = sigma2_star*length(cols_re)/nvar

# remove the smallest and largest values of sigma when plotting so that the vertical scale is not as stretched
sigma.vals = sigma.vals[3:11]
ratio = ratio[3:11,]
ratio.lo = ratio.lo[3:11,]
ratio.hi = ratio.hi[3:11,]

p3 = qplot(sigma.vals, log(ratio$ridge), size=I(0.6), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_vline(xintercept=sigma2_star, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$ridge), ymax=log(ratio.hi$ridge)), width=0.002) +
  ggtitle("Ridge Regression") +
  theme(text = element_text(size=10))


p2 = qplot(sigma.vals, log(ratio$lasso), size=I(0.6), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  ggtitle("LASSO") +
  geom_errorbar(aes(ymin=log(ratio.lo$lasso), ymax=log(ratio.hi$lasso)), width=0.002) +
  theme(text = element_text(size=10))


p4 = qplot(sigma.vals, log(ratio$nn), size=I(0.6), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  ggtitle("Neural Network") +
  geom_errorbar(aes(ymin=log(ratio.lo$nn), ymax=log(ratio.hi$nn)), width=0.002) +
  theme(text = element_text(size=10))



p5 = qplot(sigma.vals, log(ratio$rf), size=I(0.6), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  ggtitle("Random Forest") +
  geom_errorbar(aes(ymin=log(ratio.lo$rf), ymax=log(ratio.hi$rf)), width=0.002) +
  theme(text = element_text(size=10))


png('plots_thm3.png', width=600, height=800, res=100)
grid.arrange(p3, p2, p4, p5, ncol=1)
dev.off()

