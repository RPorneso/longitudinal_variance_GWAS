simtest.RIShet2 <- function(N, t, maf, ranfx, betas, expdisp, nsim){ # CHANGED disps to expdisp
  
  # Generate and save data with file name stem dgm1e_ in local drive
  
  dat <- loglindisp2(N, t, maf, ranfx, betas, expdisp) # CHANGED FROM nonlindisp2 AND disps to expdisp
  filename <- paste0("/Users/ralphp/trajsimdat/dgm1e_", N, t, "_(", ranfx[[1]],",",ranfx[[2]],",",ranfx[[3]], ")_", "(", expdisp, ")_", nsim, ".csv") # ADDED e in filename                      
  fwrite(dat, filename)    
  
  # Extract simulated parameters of interest to compare results against (1 below is for resid at L1)
  # Now only fits RIShet model following discussion with Eivind and Espen
  
  simulated <- c(betas[[1]],betas[[2]],betas[[3]],1,ranfx[[1]],ranfx[[2]],expdisp)
  
  # Test against DGM
  # Important Note! This can only accommodate up to N = 100 since NLME only can model hetsced ind in varIdent. It is not possible to use varExp here.
  
  mod <- tryCatch(lme(phe ~ time + snp, data = dat, random = ~1|id, weights = varExp(form = ~snp), method = "ML", control = lmeControl(opt = "optim")), error = function (e) NULL)
  results <- tryCatch(c(mod$coefficients$fixed[[1]], mod$coefficients$fixed[[2]], mod$coefficients$fixed[[3]], mod$sigma, as.numeric(VarCorr(mod)[3]), NA, mod$modelStruct$varStruct[[1]]), error = function (e) NULL)
  simulated - results
  
}

# , "(", disps[[1]],",",disps[[2]], ")_" # ADD THIS BACK WHEN nonlindisp2 is used
# disps[[1]]+1,disps[[2]]+1 # ADD THIS BACK WHEN nonlindisp2 is used
# varIdent(value = c("1" = 1, "2" = 1), form = ~1|snp)  # ADD THIS BACK WHEN nonlindispersion2 IS USED AND REPLACE weights in mod line
# + 1, mod$modelStruct$varStruct[[2]]  + 1 # ADD THIS BACK WHEN nonlindisp2 is used