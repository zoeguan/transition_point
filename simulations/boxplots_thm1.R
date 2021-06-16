library(gridExtra)
library(ggplot2)
library(reshape2)

load('transition_point_thm1_lme.RData')
results.lme = results



results.file = "transition_point_thm1.RData"
load(results.file)
source('transition_point_fns.R')
source('sim_fns.R')

load("hyperparams_thm1.RData")
sigma2_bar = tau_range(edat_train_misspec, edat_test_misspec, f_train, f_test, Z_train, lambda=0, lambdak=rep(0, length(edat_train_misspec)), sigma_eps=sigma_eps)
sigma2_star.misspec = (ncol(edat_train_misspec[[1]])/ncol(Z_train[[1]]))*sigma2_bar

for (j in 1:length(sigma.vals)) {
  results[[j]][, c("lme", "lme2", "meta.re")] = results.lme[[j]][, c("lme", "lme2", "meta.re")]

}

### MSE ratio plots

load("hyperparams_thm1.RData")

### MSE boxplots 
# ratio of MSE of each model to MSE of true model for specific values of sigma

for (j in 1:length(sigma.vals)) {
  r = results[[j]]
  res = colMeans(r, na.rm=T)
  err[j, ] = res
}
err$sigma2 = sigma.vals

r1 = results[[1]]/results[[1]]$lme2
names(r1) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA", "MA2")
melted1 = melt(r1[, c(3:4, 7:8, 5:6, 9:12)])
names(melted1) = c("Model", "MSPE")
melted1$log.MSPE = log(melted1$MSPE)

r2 = results[[5]]/results[[5]]$lme2
names(r2) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA", "MA2")
melted2 = melt(r2[, c(3:4, 7:8, 5:6, 9:12)])
names(melted2) = c("Model", "MSPE")
melted2$log.MSPE = log(melted2$MSPE)

r3 = results[[13]]/results[[13]]$lme2
names(r3) = c("LME", "LME2", "LS,M", "LS,E", "L,M", "L,E", "R,M", "R,E", "RF,M", "RF,E", "NN,M", "NN,E", "MA", "MA2")
melted3 = melt(r3[, c(3:4, 7:8, 5:6, 9:12)])
names(melted3) = c("Model", "MSPE")
melted3$log.MSPE = log(melted3$MSPE)


conversion.ratio = ncol(Z_train[[1]])/ncol(edat_train_misspec[[1]])

sigma.vals = sigma.vals*conversion.ratio

g1 <- ggplot(melted1, aes(x = Model, y = log.MSPE)) +
  geom_boxplot() + ggtitle(expression(paste(bar(sigma^2), "=0"))) +
  ylab(expression(paste("log", bgroup("(", frac( "MSPE"["Model"], "MSPE"["LME"]), ")" ) )))
g1

g2.title <- paste0("=", round(sigma.vals[5], 2), " (ridge regression transition point)")
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

