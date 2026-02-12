loglindisp2 <- function(N, t, maf, ranfx, betas, expdisp){
  # Creates a dataframe based on the DGM documented in RIShet(loglindisp).Rmd
  # @Params: N (numeric) is the number of individuals
  # @Params: t (numeric) is the number of time points
  # @Params: maf (numeric) is the minor allele frequency
  # @Params: ranfx (numeric vector) contains the var of u0, u1 and their cov
  # @Params: betas (numeric vector) effect of intercept, time and snp
  # @Params: expdisp (integer) is the hetsced effect of snp in loglin var
  
  id <- rep(1:N, times = t)
  snp <- rep(rbinom(N, 2, maf), times = t)
  time <- rep(1:t, each = N)
  Omega = matrix(c(ranfx[1]^2, ranfx[3], ranfx[3], ranfx[2]^2), 2, 2)
  print(Omega)
  raneffs <- mvrnorm(N, c(0,0), Omega)
  sigma2_ij <- exp(0 + expdisp * snp) # changed per espen, but now cant fit full/RIShet model
  u0_j <- raneffs[id,1]
  u1_j <- raneffs[id,2]
  a <- (betas[1] + u0_j) + ((betas[2] + u1_j) * time) + (betas[3] * snp )
  Y <- rnorm(N*t, a, sigma2_ij)
  data.frame(phe = Y, id = id, time = time, snp = snp)
}

