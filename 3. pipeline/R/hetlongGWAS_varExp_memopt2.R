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
inputfile <- args[4]
load(inputfile)
# df <- eval(parse(text = ls()))
# df$time <- ifelse(df$LEVEL=="5", 0, ifelse(df$LEVEL=="8", 3,4))
rm(list=setdiff(ls(),"df"))
print(dim(df))

# Retrieve args from bash then calculate index & load plink files

args = commandArgs(trailingOnly = TRUE)
nblocks = as.integer(args[1])
size = as.integer(args[2])
maxcol = as.integer(args[3])
taskid = Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid = as.integer(taskid)

source("/tsd/p805/data/durable/projects/ralphp/AP/code/block_index.R")
index = block_index(taskid, nblocks, size, maxcol)

genofile <- args[5]
out <- args[6]

plink <- snp_attach(genofile)
print(dim(plink$genotypes))
print(dim(plink$fam))
print(dim(plink$map))

# Subset genotypes using index

begin = index[1]
last = index[2]

# For each SNP, subset .bed to those with phenotypes, calculate MAF, repeat variant t times, and remove NAs while creating temp df, dfa; run NLME + LRT on dfa

t1 <- Sys.time()

runmods <- function (i) {
  
  SNP <- plink$genotypes[,i]
  SNP <- subset(SNP, plink$fam$sample.ID %in% df$IID)
  
  nmaf <- length(SNP) - sum(is.na(SNP))
  smaf <- sum(SNP, na.rm = TRUE)
  maf <- smaf/(nmaf * 2)
  
  repSNP <- rep(SNP, each = 3)
  indna <- is.na(repSNP)
  dfa <- df[!indna,]
  dfa$SNP <- repSNP[!indna]
  
  mod0 <- lme(resZ_SCORE ~ time + Sex + genotyping_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dfa, random = list(IID = ~1), method = "ML", control = lmeControl(opt="optim")) # base model
  mod <- lme(resZ_SCORE ~ time + Sex + genotyping_batch_num + SNP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dfa, random = list(IID = ~1), weights = varExp(form = ~ SNP), method = "ML", control = lmeControl(opt="optim"))
  lrt <- anova(mod0, mod)
  
  dispse <- tryCatch(intervals(mod, which = "var-cov")$varStruct[[2]], error=function(e) 0)
  if(dispse != 0){dispse <- (intervals(mod, which = "var-cov")$varStruct[[2]] - intervals(mod, which = "var-cov")$varStruct[[1]]) / 1.96}
 
  c(round(maf,2), summary(mod)$tTable[[5]], summary(mod)$tTable[[20]], summary(mod)$tTable[[65]], lrt[[9]][[2]], nrow(dfa), summary(mod)$modelStruct$varStruct[[1]], dispse) # includes SE
}

estimates <- list()
estimates <- lapply(begin:last, function(i) tryCatch(runmods(i), error=function(e) NULL))

nr <- max(unlist(lapply(1:length(estimates), function(i) length(estimates[[i]]))))
estimates <- lapply(1:length(estimates), function(i) if(is.null(estimates[[i]])){rep(0,nr)} else(estimates[[i]])) # replace NULL with 0s so unlist doesn't mess up order when transposed

t2 <- Sys.time()

print(t2-t1)

# Create results table using plink$map

t3 <- Sys.time()

results <- plink$map[begin:last,]

temp <- t(matrix(unlist(estimates), nr, length(estimates))) # row (coefs) x col (M) transposed becomes M x coef
colnames(temp) <- c("maf", "beta", "se", "p", "lrtp", "n", "disp", "dispse") # includes SE for disp eff
results <- cbind(results, temp)

t4 <- Sys.time()

print(t4-t3)

# Write results in a text file

t5 <- Sys.time()

fwrite(results, paste0(out, last, ".txt"))

t6 <- Sys.time()

print(t6-t5)

quit(save="no")






















