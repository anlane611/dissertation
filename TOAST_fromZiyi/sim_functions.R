library(matrixStats)
library(sirt)
library(gtools)

## get beta distribution parameters
betaMOM = function(mm, vv){
  tmp = mm*(1-mm) / vv - 1
  tmp[tmp<0] = 0.01
  alpha = mm * tmp
  beta = (1-mm) * tmp
  return(cbind(alpha, beta))
}

## get methylation profiles of pure cell types
getOnePureRefPanel = function(pure_base, pure_sd){
  N_feature = dim(pure_base)[1] # number of features
  L = dim(pure_base)[2] # number of pure tissues
  
  tissue = matrix(0,N_feature,L)
  for(i in 1:L){
    param = betaMOM(pure_base[,i],pure_sd[,i]^2)
    tissue[,i] = rbeta(N_feature,param[,1],param[,2])
  }
  tissue[tissue < 0] = 0.01
  tissue[tissue > 1] = 0.99
  
  return(tissue)
}

## get proportion matrix 
getProportion <- function(N_sample, L = 4){
  prop.matrix.ctr = matrix(runif(N_sample*L,0.05,0.95),nrow = N_sample,ncol = L)
  prop.matrix.ctr = t(apply(prop.matrix.ctr,1,function(x) x/sum(x)))
  return(prop.matrix.ctr)
}

getProportion <- function(N_sample, cc = 100){
  require("gtools")

  # alpha.ctr=c(0.97, 4.71, 0.50, 0.35)
  alpha.ctr = c(0.1, 0.2, 0.3, 0.4)
  
  prop.matrix.ctr = rdirichlet(N_sample, alpha.ctr*cc)
  
  return(prop.matrix.ctr)
}

## get mixture
getSampleMix <- function(N_sample, pure_base, pure_sd, noise_sd = 0.1){
  K = ncol(pure_base) ## number of cell types
  p = nrow(pure_base) ## number of CpG sites
  
  ## get proportions
  trueProp = getProportion(N_sample, L = K)
  alltmp = matrix(0, p, N_sample*K)
  
  ## get mix
  obs.Y = matrix(0, p, N_sample)
  
  for(n in 1:N_sample){
    tmp = getOnePureRefPanel(pure_base, pure_sd)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, noise_sd)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }
  obs.Y[obs.Y < 0] = 0.01
  obs.Y[obs.Y > 1] = 0.99
  rownames(obs.Y) = rownames(alltmp) = rownames(tmp) = rownames(pure_base)
  return(list(obs.Y = obs.Y, trueProp = trueProp,trueMethy = alltmp,trueMethy1 = tmp))
}

