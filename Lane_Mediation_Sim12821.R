rm(list=ls())

options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

if(length(args)==0){
  print("No arguments supplied")

  n_sim <- 10          #number of replicates in this script
  i_seed <- 1          #seed that changes for each replicate for new observed data
  myseed <- 123        #seed for selecting CpG sites, mediating sites, mediating cell types
  
  Nsample <- 500       #total sample size
  Ncell <- 4           #number of cell types
  Ncpg <- 10           #number of CpG sites
  
  med_num_M <- 2       #number of CpG sites that are mediators
  med_num_cell <- 2    #number of cell types that are mediators
  med_pi <- 4          #which cell is mediator
  
  is_pi_med <- FALSE   #designates whether or not pi is a mediator
  is_M_med <- TRUE     #designates whether or not M is a mediator

  alpha <- 0.00001     #significance threshold for beta/TOAST to move into EM
  
  sigma <- 1         #noise for M
  gammadenom <- 2      #noise for Y
  

} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#-------------------------
# Setup parameters
#-------------------------


library(JuliaCall)
julia <- JuliaCall::julia_setup("/apps/julia-1.5.0/bin")
#julia <- JuliaCall::julia_setup("/Applications/JuliaPro-1.5.2-1.app/Contents/Resources/julia/Contents/Resources/julia/bin")
julia_library("Distributed")

# Here we create our parallel julia processes
julia_command("addprocs(7)")

julia_command("@everywhere using Distributions")
julia_command("@everywhere using LinearAlgebra")
julia_command("@everywhere using SharedArrays")
julia_command("@everywhere using Random")

library(RcppEigen)
library(TOAST)
#source("TOAST_fromZiyi/fitModel.R")
source("fitModel.R")
library(RefFreeEWAS)
library(EpiDISH)
library(TCA)
library(boot)
library(MICS)
source("mediation_sim_functions_1120.R")
#load("TOAST_fromZiyi/BloodPureDNAmethy_4ct.rda")
load("BloodPureDNAmethy_4ct.rda")








# n_sim <- 2          #number of replicates in this script
# i_seed <- 1          #seed that changes for each replicate for new observed data
# myseed <- 123        #seed for selecting CpG sites, mediating sites, mediating cell types
# 
# Nsample <- 500       #total sample size
# Ncell <- 2           #number of cell types
# Ncpg <- 10           #number of CpG sites
# 
# med_num_M <- 1       #number of CpG sites that are mediators
# med_num_cell <- 1    #number of cell types that are mediators
# med_pi <- 2          #which cell type is the mediator (for testing purposes 2/8)
# 
# is_pi_med <- FALSE   #designates whether or not pi is a mediator
# is_M_med <- TRUE     #designates whether or not M is a mediator
# 
# alpha <- 0.000001    #significance threshold for beta/TOAST to move into EM
# 
# sigma <- 0.01         #noise for M
# gammadenom <- 5      #noise for Y




pure.mean0 = pure.mean[,1:Ncell]
pure.sd0 = pure.sd[,1:Ncell]

#--------------------------------------
# generate exposure and mediator sites
#--------------------------------------

#exposure vector
set.seed(myseed)
covar1 <- matrix(rnorm(Nsample),nrow=Nsample,ncol=1)
Ncov <- ncol(covar1)
probs <- plogis(3*covar1)
E <- rbinom(Nsample,1,probs)
#E <- c(rep(1,Nsample/2),rep(0,Nsample/2))
E.exp <- which(E==1)
E.unexp <- which(E==0)


#set which cell types are mediators
#med.pi <- sort(sample(seq(1,Ncell),med_num_cell))
med.pi <- med_pi
nonmed.pi <- setdiff(seq(1,Ncell),med.pi)

sig.list <- comp.sig <- vector("list", n_sim)
EMpvals_hybrid <- TCApvals <- pval_MultiMed <- matrix(NA,n_sim,Ncell)

beta <- matrix(NA,n_sim,Ncell)

tau <- sigma

 for (sim in 1:n_sim) {
#   #
   print(sim)

  set.seed(i_seed*1000+sim)

  rows <- which(abs(rowMeans(pure.mean0)-0.5) < 0.1)
  rows <- sample(rows, Ncpg)
  pure.mean <- pure.mean0[rows,]
  pure.sd = pure.sd0[rows,]
  pure.sd[pure.sd == 0] = 1e-6 ## avoid NA for estimating parameters of beta distribution
  
  
  #set which CpG sites are mediators
  #med.sites <- sort(sample(seq(1,nrow(pure.mean)),round(med_perc_M*nrow(pure.mean))))
  #med.sites <- sort(sample(seq(1,nrow(pure.mean)),med_num_M))
  med.sites <- sample(seq(1,nrow(pure.mean)),1)
  nonmed.sites <- setdiff(seq(1,nrow(pure.mean)),med.sites)
  
  ## Generate mediation effect size vector for sites, this links trueMethy to E.
  ## I'll make one non-zero entry for each site
  beta_effect <- runif(1,0.1,0.5)
  med.exp.M <- matrix(0, nrow=length(med.sites), ncol=Ncell)
  for(i in 1:nrow(med.exp.M)) {
    med.exp.M[i, med.pi] = beta_effect
  }
  beta[n_sim,] <- med.exp.M
  
  ## changed 2/8 so that multiple sites can have same cell type mediator
  theta_effect <- runif(1,0.1,0.5)
  cpg.celltype.thetas <- matrix(0,length(med.sites),Ncell)
  for(i in 1:nrow(cpg.celltype.thetas)) {
    cpg.celltype.thetas[i, med.pi] = theta_effect
  }
  
  direct_effect <- runif(1,0.01,0.05)
  beta_cov <- runif(1,0.01,0.05)
  
  #-------------------------------------
  #  Generate Observed Data
  #-------------------------------------
  
  
  
  out = getSampleMix(N_sample = Nsample,
                     pure.mean[,1:Ncell], pure.sd[,1:Ncell],
                     sigma, tauvec,
                     E.exp, E.unexp, med.pi, nonmed.pi, med_exp_pi,
                     med.sites, nonmed.sites, med.exp.M,
                     is_pi_med, is_M_med, beta_cov, covar1)
  Y.raw = out$obs.Y ## observed DNA methylation samples
  trueProp = out$trueProp ## proportions in samples
  trueMethy = out$trueMethy ## full true Methylation (M*pi)

  #-------------------------------------
  #  Generate Outcome
  #-------------------------------------
  
  trueMethy.med <- matrix(trueMethy[med.sites,],ncol=Ncell*Nsample)
  O3 = rep(NA, Nsample)
  for(i in 1:Nsample){
    trueMethy.ind = trueMethy.med[,seq(i, i+Nsample*(Ncell-1), by=Nsample),drop=FALSE]
    tmp = trueMethy.ind
    O3[i] = sum(tmp * cpg.celltype.thetas)
  }
  
  gamma.sd <- sd(O3)/gammadenom
  
  O = O3 + direct_effect*E + rnorm(Nsample, sd=gamma.sd)
  
  
  #---------------------------------------------
  #  Run EM, bootstrap
  #---------------------------------------------
  

  
  #sites <- which(res$joint[,2]<alpha)
  #print(sites)
  #nsites <- length(sites)
  M <- matrix(Y.raw[med.sites,],ncol=Nsample)
  coefs <- getTOAST(E,O,M,trueProp,covar1)$allbetas #get initial values for beta from TOAST
  coefs.EM <- matrix(coefs,nrow=Ncell*(Ncov+2)) #will need to change w/ covariates
  Prop <- trueProp
  
  design2 <- data.frame(O=O,E = as.factor(E),covar1)
  Design_out2 <- makeDesign(design2, Prop)
  fm2 <- fitModel(Design_out2, M)
  res <- csTest(fm2, coef=c("O"), sort = FALSE)
  pvals <- matrix(unlist(res)[1:(length(unlist(res))-3)],ncol=7,byrow = TRUE)[,6]
  theta2.toast <- ifelse(pvals>0.1,0,mean(O))
  theta2.toast1 <- ifelse(pvals>0.1,rnorm(1,mean=mean(O),sd=sd(O)),mean(O))
  
  
  #get EM/bootstrap results
  # time1 <- Sys.time()
  # initmethod <- "true"
  # myres_true <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  # EMpvals_true[sim,] <- myres_true[[2]]
  # initmethod <- "rand"
  # myres_rand <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  # EMpvals_rand[sim,] <- myres_rand[[2]]
  initmethod <- "toast"
  # myres_toast <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  # EMpvals_toast[sim,] <- myres_toast[[2]]
  myres_hybrid <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, covar1, med.pi, initmethod, theta_effect, theta2.toast1)
  
  EMpvals_hybrid[sim,] <- myres_hybrid[[2]]
  
  

  # #---------------------------------------------
  # #  Competing methods
  # #---------------------------------------------

  
  time0 <- Sys.time()
  dat <- data.frame(E=E,O=O,M=t(M),prop=trueProp,covar=covar1)
  conf_level <- 1-(.05/(Ncell*med_num_M))
  myboot2 <- function(dat,i,Ncell,med_num_M,Ncov){
    Nsample <- nrow(dat)
    .dat <- dat[i,]
    Mb <- matrix(.dat[,3:(3+med_num_M-1)],nrow=med_num_M,ncol=Nsample)
    propb <- as.matrix(.dat[,(3+med_num_M):(2+med_num_M+Ncell)],nrow=Nsample,ncol=Ncell)
    covmatb <- matrix(.dat[,(Ncell+3+med_num_M):ncol(.dat)],nrow=Nsample,ncol=Ncov)
    rownames(propb) <- colnames(Mb) <- 1:Nsample
    colnames(propb) <- 1:Ncell
    rownames(Mb) <- 1:med_num_M
    
    tca.modb <- tca(Mb,propb,verbose=FALSE)
    est.Mikb <- tensor(Mb,tca.modb,verbose=FALSE) #estimated M_ik values (list of length Ncell, each of length Nsite)
    
    getIE(.dat[,2],.dat[,1],covmatb,est.Mikb)$obsIE
  }
  test2 <- boot(dat, myboot2, R = 1000, Ncell = Ncell, med_num_M = med_num_M, Ncov = Ncov, parallel = "multicore", ncpus=8)
  #tcares <- lapply(1:Ncell, function(x2) boot.ci(boot.out = test2, conf = conf_level, type = "perc", index = x2)$percent)
  TCApvals.tmp <- matrix(NA,length(med.sites),Ncell)
  for(i in 1:Ncell){
    TCApvals.tmp[,i] <- 2*min(length(which(test2$t[,i]<0)),
                           length(which(test2$t[,i]>0)))/nrow(test2$t)
  }
  TCApvals[sim,] <- TCApvals.tmp
  print(Sys.time()-time0)



  #X.cov <- matrix(rnorm(Nsample),nrow=Nsample,ncol=1)
  out <- mics(meth_data = Y.raw, S = E, X = covar1, Y = O, cell_prop = t(trueProp))
   
  pval_MultiMed[sim,] <- out$pval_joint_MultiMed[med.sites,]
  #micSig <- ifelse(pval_MultiMed[med.sites,]<.05,1,0)


} #end sim

#filename1a <- paste0("EM_truePvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
#filename1b <- paste0("EM_randPvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
#filename1c <- paste0("EM_toastPvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename1d <- paste0("EM_hybridPvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename2 <- paste0("TCA_Pvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename3 <- paste0("MICS_Pvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)

#write.table(EMpvals_true, filename1a, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(EMpvals_rand, filename1b, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(EMpvals_toast, filename1c, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(EMpvals_hybrid, filename1d, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TCApvals, filename2, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pval_MultiMed, filename3, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

#write.table(beta,"betas",append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)