library(gridExtra)
library(ggplot2)
library(reshape2)

load("transition_point_thm1_opt.RData")
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



load("transition_point_thm1.RData")



### MSE ratio plots

bar.width = .003

p1 = qplot(sigma2.vals*length(cols_re)/10, log(ratio$csl), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_vline(xintercept=tpt_ls, color="#009900", linetype = "dotted", size=1.2) +
  geom_vline(xintercept=sigma2_bar, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$csl), ymax=log(ratio.hi$csl)), width=bar.width) +
  ggtitle("Least Squares") +
  theme(text = element_text(size=9))

p2 = qplot(sigma2.vals*length(cols_re)/10, log(ratio$ridge), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_vline(xintercept=tpt_r, color="#009900", linetype = "dotted", size=1.2) +
  geom_vline(xintercept=sigma2_bar_ridge, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$ridge), ymax=log(ratio.hi$ridge)), width=bar.width) +
  ggtitle("Ridge Regression") +
  theme(text = element_text(size=9))


png('plots_thm1_opt.png', width=800, height=800, res=150)
grid.arrange(p1, p2, ncol=1)
dev.off()
