runtest <- function(dat, type){
  # updated simulation function, basis runsim2 
  # reads files from directory and loads it as dat which is fed to runtest
  # @Params: dat is the directory and filename of df with snp, time, id & pheno
  # @Params: type refers to model to be fitted: disp1 is DGLM
  #          disp2 & 3 are varID & varExp, respectively.
  
  if (type == 'mean'){
    lm(phe~snp+snp*time,data = dat)
  }
  else if (type == 'ggi'){
    leveneTest(dat$phe~as.factor(dat$snp)*as.factor(dat$time), center = median)
  }
  else if (type == 'gei'){
    leveneTest(dat$phe~as.factor(dat$snp)*as.factor(dat$time), center = median)
  }
  else if (type == 'disp1'){
    dglm(dm_phe~time+snp, dformula = ~snp, method = "ML", zkeep = TRUE, data = dat)
  }
  else if (type == 'disp2'){ 
    base <- lme(phe ~ time, data = dat, random = list(id=~1), method = "ML", control = lmeControl(opt = "optim"))
    mod <- lme(phe ~ time + snp, data = dat, random = list(id=~1), weights = varIdent(value = c("1"=1, "2"=1), form = ~1|snp), method = "ML", control = lmeControl(opt = "optim")) # value ensures snp 0 is the ref
    lrt <- anova(base, mod) # cannot use this if there is mean effect
    lrt$'p-value'[2]
  }
  else if (type == 'disp3'){
    base <- lme(phe ~ time, data = dat, random = list(id=~1), method = "ML", control = lmeControl(opt = "optim"))
    mod <- lme(phe ~ time + snp, data = dat, random = list(id=~1), weights = varExp(form = ~ snp), method = "ML", control = lmeControl(opt = "optim"))
    lrt <- anova(base, mod) # cannot use this if there is mean effect
    lrt$'p-value'[2]
  }
}
