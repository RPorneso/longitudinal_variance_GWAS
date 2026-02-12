rm(list=ls())
library(nlme)
library(parallel)
library(MASS)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
source("~/Desktop/AP Genetics/loglindisp2.R")
source("~/Desktop/AP Genetics/loglindispersion2.R")
source("~/Desktop/AP Genetics/simtest.RIShete3.R")
source("~/Desktop/AP Genetics/simtest.RISRIhete3.R")
source("~/Desktop/AP Genetics/parsimtest.RIShete3.R")
source("~/Desktop/AP Genetics/parsimtest.RISRIhete3.R")
source("~/Desktop/AP Genetics/runtest.R")

# generate pheno based on dgmA, B, C under different exp disp effect
# loglindisp2.R, loglindispersion2.R nested in simtest.RIShete3 and simtest.RISRIhete3 nested in parsimtest.RIShete3 and parsimtest.RISRIhete3

datAe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0), betas = c(0,0,0), expdisp = 0, nsim = 100) 
datAe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0), betas = c(0,0,0), expdisp = 0.01, nsim = 100) 
datAe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0), betas = c(0,0,0), expdisp = 0.05, nsim = 100) 
datAe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) 

datBe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0.2,0), betas = c(0,0,0), expdisp = 0, nsim = 100) 
datBe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0.2,0), betas = c(0,0,0), expdisp = 0.01, nsim = 100) 
datBe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0.2,0), betas = c(0,0,0), expdisp = 0.05, nsim = 100) 
datBe <- parsimtest.RIShete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0.2,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) 

datCe <- parsimtest.RISRIhete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0.2,0,0,0), betas = c(0,0,0), expdisp = 0, nsim = 100)
datCe <- parsimtest.RISRIhete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0.2,0,0,0), betas = c(0,0,0), expdisp = 0.01, nsim = 100) 
datCe <- parsimtest.RISRIhete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0.2,0,0,0), betas = c(0,0,0), expdisp = 0.05, nsim = 100) 
datCe <- parsimtest.RISRIhete3(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0.2,0,0,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) 

# TEST USING LOGLIN MODEL VIA LRT AND RETURN PVALS
# runtest.R (copy approach in powsim-RIhet(e).R)

dats <- list.files("/Users/ralphp/trajsimdat/power/loglin", pattern = ".csv", full.names = TRUE)
dfs <- mclapply(dats, fread, mc.cores = 7)

# calculate lrt pvalue from comparing base and mod to detect hetsced effect

dats[c(1,101,201,301,401,501, 601, 701, 801, 901, 1001, 1101)] # check order to create labels (next 2 lines)
dgm <- rep(c("RIRIe", "RIe", "RISe"), each = 100*4) # 100 is nsim and 4 is for the 4 exp disp effects
tau <- rep(c(0.01,0.05,0.1,0), each = 100) # 100 is nsim 

pvals_e <- mclapply(1:length(dfs), function(x)
{
  dat <- dfs[[x]] 
  runtest(dat,type="disp3") # varExp
}, mc.cores = 7)

save(pvals_e, file = "/Users/ralphp/trajsimdat/power/pvals_e.Rda") # moved it to AP Genetics manually

p <- unlist(pvals_e)
res_e <- cbind.data.frame(dgm,tau,p)

save(res_e, file = "/Users/ralphp/trajsimdat/power/res_e.Rda") # moved it to AP Genetics manually

power <-  res_e %>%
  mutate(success = if_else(p < .05, 1, 0)) %>%
  group_by(dgm, tau) %>%
  summarise(pow = sum(success) / 100, .groups = 'drop') # 100 is nsim

# PLOT

ggplot(power, aes(x = tau, y = pow, color = dgm, group = dgm)) + 
  geom_line(linetype = 1) +
  geom_point() +
  geom_hline(yintercept = .80, linetype = 3) +
  geom_hline(yintercept = .05, linetype = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title=element_text(size=20, h = 0.5)) +
  labs(x = expression(tau), y = "Power") +
  ggtitle("longhetGWAS_loglin (N=10000, t=3)")
