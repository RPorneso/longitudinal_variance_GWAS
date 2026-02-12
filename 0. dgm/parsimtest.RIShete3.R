parsimtest.RIShete3 <- function(N, t, maf, ranfx, betas, expdisp, nsim) {
  mclapply(1:nsim, FUN = function(x) simtest.RIShete3(N, t, maf, ranfx, betas, expdisp, x), mc.cores = 7) # CHANGED disps to expdisp
}