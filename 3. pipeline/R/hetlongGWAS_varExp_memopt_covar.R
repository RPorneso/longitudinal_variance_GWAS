#!/usr/bin/Rscript

rm(list=ls())

local({
  r = getOption("repos")
  r["CRAN"] <- "https://cran.tsd.usit.no"
  r["BIOCONDUCTOR"] <- "https://bioconductor.tsd.usit.no"
  r["BIOCONDUCTOR-annotation"] <- "https://bioconductor.tsd.usit.no/data/annotation"
  r["BIOCONDUCTOR-experiment"] <- "https://bioconductor.tsd.usit.no/data/experiment"
  r["BioC_mirror"] <- "https://bioconductor.tsd.usit.no"
  options(repos=r)
  options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)
  options(download.file.method = "libcurl")
  options("install.lock"=FALSE)
  options(bitmapType="cairo")
})

LIB = "/cluster/projects/p805/ralph/software/R"
.libPaths(c(.libPaths(),LIB))

installed_packages <- rownames(installed.packages())

if(!( "MASS" %in% installed_packages)){
  install.packages(
    "https://cran.tsd.usit.no/src/contrib/Archive/MASS/MASS_7.3-60.tar.gz",
    repos=NULL,type="source",lib=LIB
  )
}

if(!( "Matrix" %in% installed_packages)){
  install.packages(
    "https://cran.tsd.usit.no/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz",
    repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB
  )
}
library(Matrix)

if(!( "bigstatsr" %in% installed_packages)){
  install.packages(
    # "https://cran.tsd.usit.no/src/contrib/Archive/bigstatsr/bigstatsr_1.5.6.tar.gz"
    "/cluster/projects/p805/crayner/software/bigstatsr_1.5.12.tar.gz",
    repos=NULL,type="source",lib=LIB
  )
}
library(bigstatsr)

if(!( "bigsnpr" %in% installed_packages)){install.packages("https://cran.tsd.usit.no/src/contrib/Archive/bigsnpr/bigsnpr_1.11.4.tar.gz",repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB)}
library(bigsnpr)

if(!( "data.table" %in% installed_packages)){install.packages("https://cran.tsd.usit.no/src/contrib/Archive/data.table/data.table_1.15.2.tar.gz",repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB)}
library(data.table)

if(!( "dplyr" %in% installed_packages)){install.packages("https://cran.tsd.usit.no/src/contrib/Archive/dplyr/dplyr_1.1.3.tar.gz",repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB)}
library(dplyr)

if(!( "tidyr" %in% installed_packages)){install.packages("https://cran.tsd.usit.no/src/contrib/Archive/tidyr/tidyr_1.3.0.tar.gz",repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB)}
library(tidyr)

if(!( "nlme" %in% installed_packages)){install.packages("https://cran.tsd.usit.no/src/contrib/Archive/nlme/nlme_3.1-166.tar.gz",repos=NULL,type="source",verbose=T,dependencies=T, lib=LIB)}
library(nlme)

library(parallel) # part of base R, no need to install from cran.tsd

# Retrieve args from bash & load analysis df

args = commandArgs(trailingOnly = TRUE)
inputfile <- args[1]
load(inputfile)
#df <- eval(parse(text = ls()))
#df$time <- ifelse(df$LEVEL=="5", 0, ifelse(df$LEVEL=="8", 3,4))
rm(list=setdiff(ls(),"df"))
print(dim(df))

# Retrieve args from bash then calculate index & load plink files

args = commandArgs(trailingOnly = TRUE)
taskid = Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid = as.integer(taskid)
i = taskid

genofile <- args[2]
snpindex <- args[3]
out <- args[4]

plink <- snp_attach(genofile)
print(dim(plink$genotypes))
print(dim(plink$fam))
print(dim(plink$map))

sigsnps <- as.vector(fread(snpindex))
print(sigsnps)

# For each SNP, subset .bed to those with phenotypes, calculate MAF, repeat variant t times, and remove NAs while creating temp df, dfa; run NLME + LRT on dfa

t1 <- Sys.time()

runmods <- function (i) {
  
  SNP <- plink$genotypes[,sigsnps[[i]]]
  SNP <- subset(SNP, plink$fam$sample.ID %in% df$IID)
  
  nmaf <- length(SNP) - sum(is.na(SNP))
  smaf <- sum(SNP, na.rm = TRUE)
  maf <- smaf/(nmaf * 2)
  
  repSNP <- rep(SNP, each = 3)
  indna <- is.na(repSNP)
  dfa <- df[!indna,]
  dfa$SNP <- repSNP[!indna]
  
  mod0 <- lme(resZ_SCORE ~  time + Sex + genotyping_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dfa, random = list(IID = ~1), weights = varComb(varExp(form = ~ PC1), varExp(form = ~ PC2), varExp(form = ~ PC3), varExp(form = ~ PC4), varExp(form = ~ PC5), varExp(form = ~ PC6), varExp(form = ~ PC7), varExp(form = ~ PC8), varExp(form = ~ PC9), varExp(form = ~ PC10)), method = "ML", control = lmeControl(opt="optim")) # base model SHOULD INCLUDE AARGANG FOR READ/MATH BUT NOT LOGBMI/HEIGHT
  mod <- lme(resZ_SCORE ~ time + Sex + genotyping_batch_num + SNP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dfa, random = list(IID = ~1), weights = varComb(varExp(form = ~ SNP), varExp(form = ~ PC1), varExp(form = ~ PC2), varExp(form = ~ PC3), varExp(form = ~ PC4), varExp(form = ~ PC5), varExp(form = ~ PC6), varExp(form = ~ PC7), varExp(form = ~ PC8), varExp(form = ~ PC9), varExp(form = ~ PC10)), method = "ML", control = lmeControl(opt="optim")) # SHOULD INCLUDE AARGANG FOR READ/MATH BUT NOT LOGBMI/HEIGHT
  lrt <- anova(mod0, mod)
  
  c(round(maf,2), summary(mod)$tTable[[5]], summary(mod)$tTable[[20]], summary(mod)$tTable[[65]], lrt[[9]][[2]], nrow(dfa), summary(mod)$modelStruct$varStruct[[1]]) # CHANGE 5, 20, 65 TO 6, 22, 70 FOR READ/MATH
}

estimates <- list()
estimates <- lapply(i, function(i) tryCatch(runmods(i), error=function(e) NULL))

t2 <- Sys.time()

print(t2-t1)

# Create results table using plink$map

t3 <- Sys.time()

results <- plink$map[sigsnps[[i]],]

temp <- t(matrix(unlist(estimates), length(estimates[[1]]), length(estimates))) # row (coefs) x col (M) transposed becomes M x coef
colnames(temp) <- c("maf", "beta", "se", "p", "lrtp", "n", "disp")
results <- cbind(results, temp)

t4 <- Sys.time()

print(t4-t3)

# Write results in a text file

t5 <- Sys.time()

fwrite(results, paste0(out, i, ".txt"))

t6 <- Sys.time()

print(t6-t5)

quit(save="no")






















