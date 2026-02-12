simtest.RISRIhet2 <- function(N, t, maf, ranfx, betas, expdisp, nsim){ # CHANGED disps TO expdisp
  
  # Generate and save data with file name stem dgm2e_ in local drive
  
  dat <- loglindispersion2(N, t, maf, ranfx, betas, expdisp) 
  filename <- paste0("/Users/ralphp/trajsimdat/dgm2e_", N, t, "_(", ranfx[[1]],",",ranfx[[2]],",",ranfx[[3]],",",ranfx[[4]],",",ranfx[[5]],",",ranfx[[6]], ")_", "(", expdisp, ")_", nsim, ".csv") # ADDED e IN FILENAME                       
  fwrite(dat, filename)   
  
  # Extract simulated parameters of interest to compare results against (1 below is for resid at L1)
  
  simulated <- c(betas[[1]],betas[[2]],betas[[3]],1,ranfx[[1]],ranfx[[2]],expdisp)
  
  # Test against DGM
  # Important Note! This can only accommodate up to N = 100 since NLME only can model hetsced ind in varIdent. It is not possible to use varExp here.
  # Now only fit RIhet model following conversation with Eivind and Espen
  
  mod <- tryCatch(lme(phe ~ time + snp, data = dat, random = ~1|id, weights = varExp(form = ~snp), method = "ML", control = lmeControl(opt = "optim")), error = function (e) NULL)
  results <- tryCatch(c(mod$coefficients$fixed[[1]], mod$coefficients$fixed[[2]], mod$coefficients$fixed[[3]], mod$sigma, as.numeric(VarCorr(mod)[3]), NA, mod$modelStruct$varStruct[[1]]), error = function (e) NULL)
  RIhet <- simulated - results
  
}
