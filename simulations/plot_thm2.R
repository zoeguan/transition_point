library(gridExtra)
library(ggplot2)
library(reshape2)

load("transition_point_thm2.RData")
source('transition_point_fns.R')
source('sim_fns.R')
load("hyperparams_thm1.RData")


ratio = data.frame(sigma.vals = sigma.vals,
                   csl=err$csl / err$merged,
                   lasso=err$lasso.csl / err$lasso.m,
                   ridge=err$ridge.csl / err$ridge.m,
                   rf=err$rf.csl / err$rf.m,
                   nn=err$nn.csl / err$nn.m)
ratio.lo = ratio
ratio.hi = ratio


### bootstrap CIs for MSE ratio comparing each CSL learner to the corresponding merged learner

for (i in 1:length(sigma.vals)) {
  res.boot = boot_ci(200, results[[i]], "merged", "csl")
  ratio.lo$csl[i] = quantile(res.boot[[3]], 0.025)
  ratio.hi$csl[i] = quantile(res.boot[[3]], 0.975)
}

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

bar.width = 0.002

# scaling factor for sigma2_star, sigma2_star.ridge, etc to make the x-axis values consistent with the definition of \overline{\sigma^2} in paper
conversion.ratio = ncol(Z_train[[1]])/ncol(edat_train_misspec[[1]])

sigma.vals = sigma.vals*conversion.ratio

# scaling factor for sigma2_bar, sigma2_bar.ridge to make the x-axis values consistent with the definition of \overline{\sigma^2} in paper
conversion.ratio.2 = ncol(edat_train[[1]])/ncol(edat_train_misspec[[1]])

p1 = qplot(sigma.vals, log(ratio$csl), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_vline(xintercept=sigma2_bar[1]*conversion.ratio.2, color='red', linetype = "longdash") +
  geom_vline(xintercept=sigma2_bar[2]*conversion.ratio.2, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$csl), ymax=log(ratio.hi$csl)), width=bar.width) +
  ggtitle("Least Squares") +
  theme(text = element_text(size=9))

p3 = qplot(sigma.vals, log(ratio$ridge), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_vline(xintercept=sigma2_bar_ridge[1]*conversion.ratio.2, color='red', linetype = "longdash") +
  geom_vline(xintercept=sigma2_bar_ridge[2]*conversion.ratio.2, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$ridge), ymax=log(ratio.hi$ridge)), width=bar.width) +
  ggtitle("Ridge Regression") +
  theme(text = element_text(size=9))


p2 = qplot(sigma.vals, log(ratio$lasso), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  ggtitle("LASSO") +
  geom_errorbar(aes(ymin=log(ratio.lo$lasso), ymax=log(ratio.hi$lasso)), width=bar.width) +
  theme(text = element_text(size=9))


p4 = qplot(sigma.vals, log(ratio$nn), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  ggtitle("Neural Network") +
  geom_errorbar(aes(ymin=log(ratio.lo$nn), ymax=log(ratio.hi$nn)), width=bar.width) +
  theme(text = element_text(size=9))


p5 = qplot(sigma.vals, log(ratio$rf), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  ggtitle("Random Forest") +
  geom_errorbar(aes(ymin=log(ratio.lo$rf), ymax=log(ratio.hi$rf)), width=bar.width) +
  theme(text = element_text(size=9)) +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") 

png('plots_thm2.png', width=600, height=600, res=100)
grid.arrange(p1, p3, p2, p4, p5, ncol=2)
dev.off()

