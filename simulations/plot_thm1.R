library(gridExtra)
library(ggplot2)
library(reshape2)

load("transition_point_thm1.RData")
source('transition_point_fns.R')
source('sim_fns.R')

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

### MSE ratio plots

load("hyperparams_thm1.RData")

# transform sigma to overline(sigma)
sigma.vals = sigma.vals*length(cols_re)/nvar
sigma2_star = sigma2_star*length(cols_re)/nvar

bar.width = .0025

p1 = qplot(sigma.vals, log(ratio$csl), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) +
  geom_vline(xintercept=sigma2_star, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$csl), ymax=log(ratio.hi$csl)), width=bar.width) +
  ggtitle("Least Squares") +
  theme(text = element_text(size=9))

p3 = qplot(sigma.vals, log(ratio$ridge), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_vline(xintercept=sigma2_bar_ridge, color='red', linetype = "longdash") +
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$ridge), ymax=log(ratio.hi$ridge)), width=bar.width) +
  ggtitle("Ridge Regression") +
  theme(text = element_text(size=9))

p2 = qplot(sigma.vals, log(ratio$lasso), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$lasso), ymax=log(ratio.hi$lasso)), width=bar.width) +
  ggtitle("LASSO") +
  theme(text = element_text(size=9))

p4 = qplot(sigma.vals, log(ratio$nn), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + 
  geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$nn), ymax=log(ratio.hi$nn)), width=bar.width) +
  ggtitle("Neural Network") +
  theme(text = element_text(size=9))

p5 = qplot(sigma.vals, log(ratio$rf), size=I(0.8), xlab=expression(bar(sigma^2)), 
           ylab=expression( paste("log", bgroup("(", frac( "MSPE"["Ensemble"], "MSPE"["Merged"]), ")" ) ))) + #geom_hline(yintercept=0, color='blue', linetype = "longdash") +
  geom_errorbar(aes(ymin=log(ratio.lo$rf), ymax=log(ratio.hi$rf)), width=bar.width) +
  ggtitle("Random Forest") +
  theme(text = element_text(size=9))


png('plots_thm1.png', width=600, height=800, res=100)
grid.arrange(p1, p3, p2, p4, p5, ncol=1)
dev.off()


### MSE boxplots 
# ratio of MSE of each model to MSE of true model for specific values of sigma

for (j in 1:length(sigma.vals)) {
  r = results[[j]]
  res = colMeans(r, na.rm=T)
  err[j, ] = res
}
err$sigma2 = sigma.vals

r1 = results[[1]]/results[[1]]$lme2
names(r1) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted1 = melt(r1[, c(3:4, 7:8, 5:6, 11:12, 14)])
names(melted1) = c("Model", "MSPE")
melted1$log.MSPE = log(melted1$MSPE)

r2 = results[[5]]/results[[5]]$lme2
names(r2) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted2 = melt(r2[, c(3:4, 7:8, 5:6, 11:12, 14)])
names(melted2) = c("Model", "MSPE")
melted2$log.MSPE = log(melted2$MSPE)

r3 = results[[13]]/results[[13]]$lme2
names(r3) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted3 = melt(r3[, c(3:4, 7:8, 5:6, 11:12, 14)])
names(melted3) = c("Model", "MSPE")
melted3$log.MSPE = log(melted3$MSPE)


g1 <- ggplot(melted1, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(expression(paste(bar(sigma^2), "=0"))) +
  ylab(expression(paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
g1

g2.title <- paste0("=", round(sigma.vals[5], 2), " (least squares transition point)")
g2 <- ggplot(melted2, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(bquote(paste(bar(sigma^2)) ~ .(g2.title))) +
  ylab(expression( paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
g2

g3.title <- paste0("=", round(sigma.vals[13], 2), "")
g3 <- ggplot(melted3, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(bquote(paste(bar(sigma^2)) ~ .(g3.title))) +
  ylab(expression( paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
g3


png('boxplots_thm1.png', width=650, height=600, res=100)
grid.arrange(g1, g2, g3, ncol=1)
dev.off()


# adding random forest and univariate meta-analysis results

r1 = results[[1]]/results[[1]]$lme2
names(r1) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted1 = melt(r1[, c(3:4, 7:8, 5:6, 11:12, 9:10, 14, 13)])
names(melted1) = c("Model", "MSPE")
melted1$log.MSPE = log(melted1$MSPE)

r2 = results[[5]]/results[[5]]$lme2
names(r2) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted2 = melt(r2[, c(3:4, 7:8, 5:6, 11:12, 9:10, 14, 13)])
names(melted2) = c("Model", "MSPE")
melted2$log.MSPE = log(melted2$MSPE)

r3 = results[[13]]/results[[13]]$lme2
names(r3) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA2", "MA")
melted3 = melt(r3[, c(3:4, 7:8, 5:6, 11:12, 9:10, 14, 13)])
names(melted3) = c("Model", "MSPE")
melted3$log.MSPE = log(melted3$MSPE)


g1 <- ggplot(melted1, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(expression(paste(bar(sigma^2), "=0"))) +
  ylab(expression(paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
#g1

g2 <- ggplot(melted2, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(bquote(paste(bar(sigma^2)) ~ .(g2.title))) +
  ylab(expression( paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
#g2

g3 <- ggplot(melted3, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(bquote(paste(bar(sigma^2)) ~ .(g3.title))) +
  ylab(expression( paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
#g3


png('boxplots_thm1_rf.png', width=700, height=600, res=100)
grid.arrange(g1, g2, g3, ncol=1)
dev.off()
