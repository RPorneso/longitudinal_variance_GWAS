simtest.RISRIhete3 <- function(N, t, maf, ranfx, betas, expdisp, nsim){
  
  # Generate and save data with file name stem dgm1e_ in local drive for power simulation
  
  dat <- loglindispersion2(N, t, maf, ranfx, betas, expdisp) 
  filename <- paste0("/Users/ralphp/trajsimdat/power2/loglin/dgm1e_", N, t, "_(", ranfx[[1]],",",ranfx[[2]],",",ranfx[[3]], ")_", "(", expdisp, ")_", nsim, ".csv")                   
  fwrite(dat, filename)    
  
}

