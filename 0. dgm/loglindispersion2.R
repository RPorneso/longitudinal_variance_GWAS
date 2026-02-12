loglindispersion2 <- function(N, t, maf, ranfx, betas, expdisp){
  # Creates a dataframe based on the DGM documented in RISRIhet(loglindispersion).Rmd
  # @Params: N (numeric) is the number of individuals
  # @Params: t (numeric) is the number of time points
  # @Params: maf (numeric) is the minor allele frequency
  # @Params: ranfx (numeric vector of length 6) contains the var of u0, u1, g &
  #        : their covariances u01, u0g, u1g where g is the ranint of I in resid
  # @Params: betas (numeric vector) effect of intercept, time and snp
  # @Params: expdisp (numeric vector) is the hetsced effect of snp in loglin var
  
  id <- rep(1:N, times = t)
  snp <- rep(rbinom(N, 2, maf), times = t)
  time <- rep(1:t, each = N)
  Omega <- matrix(c(ranfx[1]^2, ranfx[4], ranfx[5], 
                    ranfx[4], ranfx[2]^2, ranfx[6], 
                    ranfx[5], ranfx[6], ranfx[3]^2), 3, 3)
  print(Omega)
  raneffs <- mvrnorm(N, c(0,0,0), Omega)
  g <- raneffs[id,3] 
  sigma2_ij <- exp(0 + expdisp * snp + g) # hetsced snp eff as loglin var
  u0_j <- raneffs[id,1]
  u1_j <- raneffs[id,2]
  a <- (betas[1] + u0_j) + ((u1_j + betas[2]) * time) + (betas[3] * snp )
  Y <- rnorm(N*t, a, sigma2_ij)
  data.frame(phe = Y, id = as.factor(id), time = time, snp = snp)
}

