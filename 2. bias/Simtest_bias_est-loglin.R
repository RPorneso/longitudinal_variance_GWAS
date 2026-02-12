rm(list = ls())
# install.packages("lqmm")
# install.packages("patchwork")
library(dplyr)
library(MASS)
library(nlme)
library(lqmm)
library(parallel)
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)
source("~/Desktop/AP Genetics/loglindisp2.R")
source("~/Desktop/AP Genetics/loglindispersion2.R")
source("~/Desktop/AP Genetics/simtest.RISRIhet2.R")
source("~/Desktop/AP Genetics/simtest.RIShet2.R")
source("~/Desktop/AP Genetics/parsimtest.RISRIhet2.R")
source("~/Desktop/AP Genetics/parsimtest.RIShet2.R")

# Simulate RIhete then RIShete, test against RIhete and return parameters of interest in lists of list - dgmAe & dgmBe (v2 with expdisp=0.1)
# Simulate RIRIhete, test against RIhete and return parameters of interest in lists of list - dgmCe (v2 with expdisp=0.1)

dgmAe_v2 <- parsimtest.RIShet2(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) # NEED TO SWITCH DIFF CALCULATION 
dgmBe_v2 <- parsimtest.RIShet2(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0.2,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) # NEED TO SWITCH DIFF CALCULATION 
dgmCe_v2 <- parsimtest.RISRIhet2(N = 10000, t = 3, maf = 0.3, ranfx = c(0.9,0,0.2,0,0,0), betas = c(0,0,0), expdisp = 0.1, nsim = 100) # NEED TO SWITCH DIFF CALCULATION 

save(dgmAe_v2, file = "~/Desktop/AP Genetics/dgmAe.Rda")
save(dgmBe_v2, file = "~/Desktop/AP Genetics/dgmBe.Rda")
save(dgmCe_v2, file = "~/Desktop/AP Genetics/dgmCe.Rda")

## Calculate mean and SE for every model fitted using each list-of-list above containing nsim sets of fixed and random effect estimates

j <- 3 # no. of dgms that completed running
models <- ls()[1:j]

xbar <- lapply(1:j, function(i) Reduce('+', get(models[[i]])) / length(get(models[[i]])))

sqdiff <- vector(mode = "list")

for (h in 1:j){
  temp <- lapply(1:length(get(models[[h]])), function(i) (get(models[[h]])[[i]] - xbar[[h]])^2)
  sqdiff[[h]] <- temp
}

SE <- lapply(1:j, function(i) sqrt(Reduce('+', sqdiff[[i]])) / length(sqdiff[[i]]))

## PLOT TOGETHER

Parameters <- c("b0", "b1", "b2", "e^alpha", "sigma^2_u0", "sigma^2_u1","tau") 
a <- cbind(Parameters, as.data.frame(xbar[[1]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED 
b <- cbind(Parameters, as.data.frame(SE[[1]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED
mydf <- cbind(a,b)
mydf <- mydf[-6,-3] # removed sigma^2_u1 which is NA
colnames(mydf) <- c("Parameters", "Bias", "stderr")

plot1 <- ggplot(mydf,aes(x=Parameters, y = Bias)) +
  geom_point(position=position_dodge(width=0.5))  +
  geom_hline(yintercept = 0, size = I(0.4), color = I("red")) +
  ylim(min = - 1.1, max = 1.1) +
  geom_errorbar(aes(x=Parameters, ymin = Bias - 1.96*stderr, ymax = Bias + 1.96*stderr),width = .1,position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text = element_text(size = 18),  axis.title.y = element_blank(), axis.text.x=element_blank(), axis.title.x = element_blank()) + 
  ggtitle("Data generation w/ (a) random intercept and (b) slope or (c) hetsced individual")

a <- cbind(Parameters, as.data.frame(xbar[[2]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED
b <- cbind(Parameters, as.data.frame(SE[[2]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED
mydf <- cbind(a,b)
mydf <- mydf[-6,-3] # removed sigma^2_u1 which is NA
colnames(mydf) <- c("Parameters", "Bias", "stderr")

plot2 <- ggplot(mydf,aes(x=Parameters, y = Bias)) +
  geom_point(position=position_dodge(width=0.5))  +
  geom_hline(yintercept = 0, size = I(0.4), color = I("red")) +
  ylim(min = - 1.1, max = 1.1) +
  geom_errorbar(aes(x=Parameters, ymin = Bias - 1.96*stderr, ymax = Bias + 1.96*stderr),width = .1,position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 18),  axis.title.y = element_text(size = 18), axis.text.x=element_blank(), axis.title.x = element_blank()) 

a <- cbind(Parameters, as.data.frame(xbar[[3]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED
b <- cbind(Parameters, as.data.frame(SE[[3]])*-1) # KEEP -1 UNTIL YOU REVERSE DIFF = SIMULATED - ESTIMATED
mydf <- cbind(a,b)
mydf <- mydf[-6,-3] # removed sigma^2_u1 which is NA
colnames(mydf) <- c("Parameters", "Bias", "stderr")

plot3 <- ggplot(mydf,aes(x=Parameters, y = Bias)) +
  geom_point(position=position_dodge(width=0.5))  +
  geom_hline(yintercept = 0, size = I(0.4), color = I("red")) +
  ylim(min = - 1.1, max = 1.1) +
  geom_errorbar(aes(x=Parameters, ymin = Bias - 1.96*stderr, ymax = Bias + 1.96*stderr),width = .1,position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 18),  axis.title.y = element_blank(), axis.title.x = element_blank()) 

plot1 / plot2 / plot3
