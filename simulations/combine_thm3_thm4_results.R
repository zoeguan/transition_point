### Theorem 4

load(paste0("thm3_thm4_RData/transition_point_thm4_", 1, "_", 1, ".RData"))
err.all = err
results.all = results

for (i in 2:13) {
  load(paste0("thm3_thm4_RData/transition_point_thm4_", i, "_", i, ".RData"))
  err.all[i, ] = err[i, ]
  results.all[[i]] = results[[i]]
}

err = err.all
results = results.all

save(err, results, scale.vals, sigma2_star, sigma2_bar, 
     scale.vals.2, cols_re,
     sigma2_star_lo, sigma2_star_hi, 
     file="transition_point_thm4.RData")



### Theorem 3

load(paste0("thm3_thm4_RData/transition_point_thm3_", 1, "_", 1, ".RData"))
err.all = err
results.all = results

for (i in 2:13) {
  load(paste0("thm3_thm4_RData/transition_point_thm3_", i, "_", i, ".RData"))
  err.all[i, ] = err[i, ]
  results.all[[i]] = results[[i]]
}

err = err.all
results = results.all

save(err, results, 
     sigma.vals, sigma2_star, cols_re,
     sigma2_bar,
     file=paste0("transition_point_thm3.RData"))


