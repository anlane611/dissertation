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
getOnePureRefPanel = function(pure_base, pure_sd, med.pi, nonmed.pi,
                              med.exp.M, med.sites, nonmed.sites, med = TRUE){
  N_feature = dim(pure_base)[1] # number of features
  L = dim(pure_base)[2] # number of pure tissues

  tissue = matrix(0,N_feature,L)
  N_med_feature = length(med.sites)
  N_nonmed_feature = length(nonmed.sites)

  ## first generate all data, without considering mediation
  for(i in 1:L) {
      param = betaMOM(pure_base[,i],pure_sd[,i]^2)
      tissue[,i] = rbeta(N_feature,param[,1],param[,2])
  }

  if (med){
      ## multiply med.exp.M for mediation sites
      tissue[med.sites,] = tissue[med.sites,] * med.exp.M

  }

  tissue[tissue < 0] = 0.01
  tissue[tissue > 1] = 0.99

  return(tissue)
}


getProportion <- function(N_sample, cc = 100, E.exp, E.unexp, med.pi,
                          nonmed.pi, med_exp_pi, med = TRUE){
  require("gtools")

  # alpha.ctr=c(0.97, 4.71, 0.50, 0.35)
  alpha.ctr = c(0.1, 0.2, 0.3, 0.4)

  prop.matrix.ctr = matrix(0,N_sample,length(alpha.ctr))
    if (med){
      prop.matrix.ctr[E.exp,] = rdirichlet(length(E.exp), alpha.ctr*cc)
      prop.matrix.ctr[E.exp,med.pi] = prop.matrix.ctr[E.exp,med.pi]*med_exp_pi
      prop.matrix.ctr[E.exp,nonmed.pi] = prop.matrix.ctr[E.exp,nonmed.pi] +
          (1-rowSums(prop.matrix.ctr)[1:length(E.exp)])/length(nonmed.pi)

      prop.matrix.ctr[E.unexp,] = rdirichlet(length(E.unexp), alpha.ctr*cc)
    } else {
      prop.matrix.ctr = rdirichlet(N_sample, alpha.ctr*cc)
  }

  return(prop.matrix.ctr)
}

## get mixture
getSampleMix <- function(N_sample, pure_base, pure_sd, noise_sd = 0.1,
                         E.exp, E.unexp, med.pi, nonmed.pi, med_exp_pi,
                         med.sites, nonmed.sites, med_exp_M,
                         is_pi_med, is_M_med){
  K = ncol(pure_base) ## number of cell types
  p = nrow(pure_base) ## number of CpG sites

  ## get proportions
  trueProp = getProportion(N_sample, cc=100, E.exp, E.unexp, med.pi,
                           nonmed.pi, med_exp_pi, is_pi_med)
  alltmp = matrix(0, p, N_sample*K)

  ## get mix
  obs.Y = matrix(0, p, N_sample)

  for(n in E.exp){ #exposed group  - relationship with mediator
    tmp = getOnePureRefPanel(pure_base, pure_sd, med.pi, nonmed.pi,
                             med.exp.M, med.sites, nonmed.sites, med = is_M_med)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, noise_sd)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }
  for(n in E.unexp){ #unexposed group - no relationship with mediator
    tmp = getOnePureRefPanel(pure_base, pure_sd, med.pi, nonmed.pi,
                             med_exp_M, med.sites, nonmed.sites, med = FALSE)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, noise_sd)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }

  obs.Y[obs.Y < 0] = 0.01
  obs.Y[obs.Y > 1] = 0.99
  rownames(obs.Y) = rownames(alltmp) = rownames(tmp) = rownames(pure_base)
  return(list(obs.Y = obs.Y, trueProp = trueProp,trueMethy = alltmp,trueMethy1 = tmp))
}

## get coefficients for k mediation models
med_coeff <- function(vec.out,E.prod){
  return(c(E.prod%*%vec.out))
}

## get indirect effect with permutation
myperm <- function(E,O,M,perm.E=TRUE,perm.O=TRUE){
  E.p <- if(perm.E){sample(E,length(E),replace = FALSE)}else{E}
  O.p <- if(perm.O){sample(O,length(O),replace = FALSE)}else{O}

  #outcome model
  out.design.p <- cbind(E.p,M)
  out.coeffs.p <- c(solve(t(out.design)%*%out.design)%*%t(out.design)%*%O.p)[2:(ncol(M)+1)]

  #mediator model
  E.prod <- solve(t(E.p)%*%E.p)%*%t(E.p)
  med.coeffs.p <- apply(M,2,function(x2) E.prod%*%x2)

  #calculate indirect effect
  perm_indirect <- abs(out.coeffs.p %*% med.coeffs.p)

  return(perm_indirect)
}
