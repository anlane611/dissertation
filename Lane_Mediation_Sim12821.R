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
  beta_effect <- 1.4   #relationship between mediator M and exposure
  theta_effect <- 0.4  #relationship between mediator M and outcome
  
  is_pi_med <- FALSE   #designates whether or not pi is a mediator
  is_M_med <- TRUE     #designates whether or not M is a mediator
  
  med_exp_pi <- 1      #relationship between mediator pi and exposure
  med_out_pi <- 1      #relationship between mediator pi and outcome
  
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

library(RcppEigen)
library(TOAST)
source("TOAST_fromZiyi/fitModel.R")
#source("fitModel.R")
library(RefFreeEWAS)
library(EpiDISH)
library(TCA)
library(boot)
library(MICS)
source("mediation_sim_functions_1120.R")
load("TOAST_fromZiyi/BloodPureDNAmethy_4ct.rda")
#load("BloodPureDNAmethy_4ct.rda")








n_sim <- 1          #number of replicates in this script
i_seed <- 1          #seed that changes for each replicate for new observed data
myseed <- 123        #seed for selecting CpG sites, mediating sites, mediating cell types

Nsample <- 500       #total sample size
Ncell <- 4           #number of cell types
Ncpg <- 10           #number of CpG sites

med_num_M <- 1       #number of CpG sites that are mediators
med_num_cell <- 1    #number of cell types that are mediators
med_pi <- 4          #which cell type is the mediator (for testing purposes 2/8)
beta_effect <- 1.2   #relationship between mediator M and exposure
theta_effect <- 0.2  #relationship between mediator M and outcome

is_pi_med <- FALSE   #designates whether or not pi is a mediator
is_M_med <- TRUE     #designates whether or not M is a mediator

med_exp_pi <- 1      #relationship between mediator pi and exposure
med_out_pi <- 1      #relationship between mediator pi and outcome

alpha <- 0.000001    #significance threshold for beta/TOAST to move into EM

sigma <- 1         #noise for M
gammadenom <- 2      #noise for Y




pure.mean0 = pure.mean[,1:Ncell]
pure.sd0 = pure.sd[,1:Ncell]

#--------------------------------------
# generate exposure and mediator sites
#--------------------------------------
set.seed(myseed)

#exposure vector
E <- c(rep(1,Nsample/2),rep(0,Nsample/2))
E.exp <- which(E==1)
E.unexp <- which(E==0)

# reduce dimension - select CpG sites close to 0.5
rows <- which(abs(rowMeans(pure.mean0)-0.5) < 0.1)
rows <- sample(rows, Ncpg)
pure.mean <- pure.mean0[rows,]

#set which CpG sites are mediators
#med.sites <- sort(sample(seq(1,nrow(pure.mean)),round(med_perc_M*nrow(pure.mean))))
#med.sites <- sort(sample(seq(1,nrow(pure.mean)),med_num_M))
med.sites <- 5
nonmed.sites <- setdiff(seq(1,nrow(pure.mean)),med.sites)


#set which cell types are mediators
#med.pi <- sort(sample(seq(1,Ncell),med_num_cell))
med.pi <- med_pi
nonmed.pi <- setdiff(seq(1,Ncell),med.pi)

## Generate mediation effect size vector for sites, this links trueMethy to E.
## I'll make one non-zero entry for each site
med.exp.M <- matrix(1, nrow=length(med.sites), ncol=Ncell)
for(i in 1:nrow(med.exp.M)) {
  med.exp.M[i, med.pi] = beta_effect
}

## changed 2/8 so that multiple sites can have same cell type mediator
cpg.celltype.thetas <- matrix(0,length(med.sites),Ncell)
for(i in 1:nrow(cpg.celltype.thetas)) {
  cpg.celltype.thetas[i, med.pi] = theta_effect
}

sig.list <- comp.sig <- vector("list", n_sim)
EMpvals_true <- EMpvals_rand <- EMpvals_toast <- matrix(NA,n_sim,Ncell)
tau <- sigma
tauvec <- rep(tau,Ncell)
pure.sd = pure.sd0[rows,]
pure.sd[pure.sd == 0] = 1e-6 ## avoid NA for estimating parameters of beta distribution



for (sim in 1:n_sim) {
  #
  print(sim)

  set.seed(i_seed*1000+sim)

  
  #-------------------------------------
  #  Generate Observed Data
  #-------------------------------------
  
  
  
  out = getSampleMix(N_sample = Nsample,
                     pure.mean[,1:Ncell], pure.sd[,1:Ncell],
                     sigma, tauvec,
                     E.exp, E.unexp, med.pi, nonmed.pi, med_exp_pi,
                     med.sites, nonmed.sites, med.exp.M,
                     is_pi_med, is_M_med)
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
  
  O = O3 + rnorm(Nsample, sd=gamma.sd)
  
  
  #---------------------------------------------
  #  Run EM, bootstrap
  #---------------------------------------------
  

  
  #sites <- which(res$joint[,2]<alpha)
  #print(sites)
  #nsites <- length(sites)
  M <- matrix(Y.raw[med.sites,],ncol=Nsample)
  coefs <- getTOAST(E,O,M,trueProp)$allbetas #get initial values for beta from TOAST
  coefs.EM <- matrix(coefs,nrow=Ncell*2) #will need to change w/ covariates
  Prop <- trueProp
  
  design2 <- data.frame(O=O,E = as.factor(E))
  Design_out2 <- makeDesign(design2, Prop)
  fm2 <- fitModel(Design_out2, M)
  res <- csTest(fm2, coef=c("O"), sort = FALSE)
  pvals <- matrix(unlist(res)[1:(length(unlist(res))-3)],ncol=7,byrow = TRUE)[,6]
  theta2.toast <- ifelse(pvals>0.1,0,mean(O))
  
  #get EM/bootstrap results
  # time1 <- Sys.time()
  initmethod <- "true"
  myres_true <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  EMpvals_true[sim,] <- myres_true[[2]]
  initmethod <- "rand"
  myres_rand <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  EMpvals_rand[sim,] <- myres_rand[[2]]
  initmethod <- "toast"
  myres_toast <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, initmethod, theta_effect, theta2.toast)
  EMpvals_toast[sim,] <- myres_toast[[2]]
  
  # print(Sys.time()-time1)
  #print(sig.list[[sim]])
  
  
  # rownames(trueProp) <- colnames(M) <- 1:Nsample
  # colnames(trueProp) <- 1:Ncell
  # rownames(M) <- 1:med_num_M
  # 
  # tca.mod <- tca(M,trueProp,verbose=FALSE)
  # est.Mik <- tensor(M,tca.mod,verbose=FALSE) #estimated M_ik values (list of length Ncell, each of length Nsite)
  # 
  # summary(fastLm(O~E+t(est.Mik[[1]])+t(est.Mik[[2]])+t(est.Mik[[3]])+t(est.Mik[[4]])))
  
  # #---------------------------------------------
  # #  Competing methods
  # #---------------------------------------------

  # dat <- data.frame(E=E,O=O,M=t(M),prop=trueProp)
  # time0 <- Sys.time()
  # myboot1 <- function(dat,i,med_num_M,Ncell){
  #   .dat <- dat[i,]
  #   getTOAST(.dat[,1],.dat[,2],t(matrix(.dat[,3:(3+med_num_M-1)])),.dat[,(3+med_num_M):(2+med_num_M+Ncell)])$indef
  #   
  # }
  # test1 <- boot(dat, myboot1, R = 1000, med_num_M = med_num_M, Ncell=Ncell, parallel = "multicore", ncpus=8)
  # toastres <- lapply(1:Ncell, function(x2) boot.ci(boot.out = test1, conf = 0.9875, type = "perc", index = x2)$percent)
  
  
  time0 <- Sys.time()
  dat <- data.frame(E=E,O=O,M=t(M),prop=trueProp)
  conf_level <- 1-(.05/(Ncell*med_num_M))
  myboot2 <- function(dat,i,Ncell,med_num_M){
    Nsample <- nrow(dat)
    .dat <- dat[i,]
    Mb <- matrix(.dat[,3:(3+med_num_M-1)],nrow=med_num_M,ncol=Nsample)
    propb <- as.matrix(.dat[,(3+med_num_M):(2+med_num_M+Ncell)],nrow=Nsample,ncol=Ncell)
    rownames(propb) <- colnames(Mb) <- 1:Nsample
    colnames(propb) <- 1:Ncell
    rownames(Mb) <- 1:med_num_M
    
    tca.modb <- tca(Mb,propb,verbose=FALSE)
    est.Mikb <- tensor(Mb,tca.modb,verbose=FALSE) #estimated M_ik values (list of length Ncell, each of length Nsite)
    
    getIE(.dat[,2],.dat[,1],est.Mikb)$obsIE
  }
  test2 <- boot(dat, myboot2, R = 1000, Ncell = Ncell, med_num_M = med_num_M, parallel = "multicore", ncpus=8)
  #tcares <- lapply(1:Ncell, function(x2) boot.ci(boot.out = test2, conf = conf_level, type = "perc", index = x2)$percent)
  TCApvals <- matrix(NA,length(med.sites),Ncell)
  for(i in 1:Ncell){
    TCApvals[,i] <- 2*min(length(which(test2$t[,i]<0)),
                           length(which(test2$t[,i]>0)))/nrow(test2$t)
  }
  print(Sys.time()-time0)




  # bootSig <- matrix(NA,nrow=1,ncol=Ncell)
  # for(i in 1:Ncell){
  #   bootSig[1,i] <- ifelse(0>tcares[[i]][4]&0<tcares[[i]][5],0,1)
  #   #bootSig[2,i] <- ifelse(0>toastres[[i]][4]&0<toastres[[i]][5],0,1)
  # }
  # 
  # ## MICS
  # 
  # #load example data
  # data(example_data)
  # 
  # #perform mics
  # out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, cell_prop = P_matr,
  #             MCP.type = "FDR", maxp_sq = TRUE)
  # 
  X.cov <- matrix(rnorm(Nsample),nrow=Nsample,ncol=1)
  out <- mics(meth_data = Y.raw, S = E, X = X.cov, Y = O, cell_prop = t(trueProp))
   
  pval_MultiMed <- out$pval_joint_MultiMed[med.sites,]
  #micSig <- ifelse(pval_MultiMed[med.sites,]<.05,1,0)
  #comp.sig[[sim]] <- rbind(bootSig,micSig)
  

} #end sim

filename1a <- paste0("EM_truePvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename1b <- paste0("EM_randPvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename1c <- paste0("EM_toastPvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename2 <- paste0("TCA_Pvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)
filename3 <- paste0("MICS_Pvals_N",Nsample,"_Ncell",Ncell,"_medpi",med.pi,"_sigma",sigma,"_tau",tau,"_gammadenom",gammadenom,"_iseed",i_seed)

write.table(EMpvals_true, filename1a, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(EMpvals_rand, filename1b, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(EMpvals_toast, filename1c, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TCApvals, filename2, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pval_MultiMed, filename3, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
