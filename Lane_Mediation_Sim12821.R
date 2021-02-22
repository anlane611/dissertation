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
  beta_effect <- 1.4   #relationship between mediator M and exposure
  theta_effect <- 0.4  #relationship between mediator M and outcome
  
  is_pi_med <- FALSE   #designates whether or not pi is a mediator
  is_M_med <- TRUE     #designates whether or not M is a mediator
  
  med_exp_pi <- 1      #relationship between mediator pi and exposure
  med_out_pi <- 1      #relationship between mediator pi and outcome
  
  alpha <- 0.00001     #significance threshold for beta/TOAST to move into EM

} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## HERE IS A TEST CHANGE FOR GITHUB!!!


#-------------------------
# Setup parameters
#-------------------------

library(JuliaCall)
#julia <- JuliaCall::julia_setup("/apps/julia-1.5.0/bin")
julia <- JuliaCall::julia_setup("/Applications/JuliaPro-1.5.2-1.app/Contents/Resources/julia/Contents/Resources/julia/bin")
julia_library("Distributed")

# Here we create our parallel julia processes
julia_command("addprocs(3)")

julia_command("@everywhere using Distributions")
julia_command("@everywhere using LinearAlgebra")
julia_command("@everywhere using SharedArrays")

library(RcppEigen)
library(TOAST)
source("fitModel.R")
library(RefFreeEWAS)
library(EpiDISH)
library(TCA)
#library(MICS)
source("mediation_sim_functions_1120.R")
load("BloodPureDNAmethy_4ct.rda")



n_sim <- 1          #number of replicates in this script
i_seed <- 1          #seed that changes for each replicate for new observed data
myseed <- 123        #seed for selecting CpG sites, mediating sites, mediating cell types

Nsample <- 500       #total sample size
Ncell <- 4           #number of cell types
Ncpg <- 10           #number of CpG sites

med_num_M <- 2       #number of CpG sites that are mediators
med_num_cell <- 1    #number of cell types that are mediators
med_pi <- 2          #which cell type is the mediator (for testing purposes 2/8)
beta_effect <- 1.4   #relationship between mediator M and exposure
theta_effect <- 0.4  #relationship between mediator M and outcome

is_pi_med <- FALSE   #designates whether or not pi is a mediator
is_M_med <- TRUE     #designates whether or not M is a mediator

med_exp_pi <- 1      #relationship between mediator pi and exposure
med_out_pi <- 1      #relationship between mediator pi and outcome

alpha <- 0.000001    #significance threshold for beta/TOAST to move into EM




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
med.sites <- sort(sample(seq(1,nrow(pure.mean)),med_num_M))
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

sig.list <- tca.sig <- tca.rank <- vector("list", n_sim)
alpha_tca <- .05/(Ncell*med_num_M)    #significance threshold for TCA models


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
                     noise_sd = 0.01,
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
  
  O = O3 + rnorm(Nsample, sd=sd(O3)/5) ## I'll keep the error small for now
  
  
  #---------------------------------------------
  #  Get TOAST estimates
  #---------------------------------------------
  
  design <- data.frame(E = as.factor(E))
  Prop <- trueProp
  colnames(Prop) <- c("CD8T","BCell","Mono","Gran")[1:Ncell] ## need colnames
  Design_out <- makeDesign(design, Prop)
  
  fm <- fitModel(Design_out, Y.raw)
  
  ## estimates for E[M|A]
  coefs <- fm$coefs ## dim: 2*Ncell X NCpG

  res <- csTest(fm, coef="E1", sort = FALSE)

  
  #---------------------------------------------
  #  Run EM, bootstrap
  #---------------------------------------------
  
  #select significant CpG sites from TOAST results
  
  #sites <- which(res$joint[,2]<alpha)
  #print(sites)
  #nsites <- length(sites)
  M <- matrix(Y.raw[med.sites,],ncol=Nsample)
  coefs.EM <- matrix(coefs[,med.sites],nrow=Ncell*2) #will need to change w/ covariates
  
  #get EM/bootstrap results
  sig.list[[sim]] <- julia_call("EMBoot", M, O, E, coefs.EM, Prop, med.pi, theta_effect)
  #print(sig.list[[sim]])
  
  # #---------------------------------------------
  # #  Competing methods
  # #---------------------------------------------

  ## Use TCA estimates in mediation models
  rownames(trueProp) <- colnames(M) <- 1:Nsample
  colnames(trueProp) <- 1:Ncell
  rownames(M) <- 1:med_num_M
  
  tca.mod <- tca(M,trueProp)
  est.Mik <- tensor(M,tca.mod)

  theta2 <- pvals <- matrix(0,med_num_M, Ncell)
  theta0 <- c()
  theta1 <- c()
  if(Ncell==4){
    for(i in 1:med_num_M){
      mod <- summary(fastLm(O~E+est.Mik[[1]][i,]+
                              est.Mik[[2]][i,]+
                              est.Mik[[3]][i,]+
                              est.Mik[[4]][i,]))
      pvals[i,] <- mod$coefficients[3:(Ncell+2),4]
      theta2[i,] <- mod$coefficients[3:(Ncell+2),1]
      theta1[i] <- mod$coefficients[2,1]
    }
  }
  
  if(Ncell==3){
    for(i in 1:med_num_M){
      mod <- summary(fastLm(O~E+est.Mik[[1]][i,]+
                              est.Mik[[2]][i,]+
                              est.Mik[[3]][i,]))
      pvals[i,] <- mod$coefficients[3:(Ncell+2),4]
      theta2[i,] <- mod$coefficients[3:(Ncell+2),1]
      theta1[i] <- mod$coefficients[2,1]
    }
  }
  
  if(Ncell==2){
    for(i in 1:med_num_M){
      mod <- summary(fastLm(O~E+est.Mik[[1]][i,]+
                              est.Mik[[2]][i,]))
      pvals[i,] <- mod$coefficients[3:(Ncell+2),4]
      theta2[i,] <- mod$coefficients[3:(Ncell+2),1]
      theta1[i] <- mod$coefficients[2,1]
    }
  }
  
  
  #identify significant sites/cell types
  tca.sig[[sim]] <- ifelse(pvals<alpha_tca,1,0)
  tca.rank[[sim]] <- matrix(rank(pvals), ncol=ncol(pvals))
  
  # ## MICS
  # 
  # #load example data
  # data(example_data)
  # 
  # #perform mics
  # #out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, cell_prop = P_matr, 
  #  #           MCP.type = "FDR", maxp_sq = TRUE)
  # 
  # out <- mics(meth_data = Y.raw, S = E, X = NULL, Y = O, cell_prop = trueProp)
  # 
  # pval_MultiMed <- out$pval_joint_MultiMed
  # 
  # #FDR threshold
  # fdr_thred <- 0.2
  # 
  # ind1 <- which(pval_MultiMed[,1] < fdr_thred)
  # 
  # ind2 <- which(pval_MultiMed[,2] < fdr_thred)
  # 
  # ind3 <- which(pval_MultiMed[,3] < fdr_thred)
  
  
  

} #end sim


total.sig.em <- Reduce("+",sig.list)
total.sig.tca <- Reduce("+",tca.sig)
avg.sig.tca <- Reduce("+",tca.rank)/length(tca.rank)

filename1 <- paste0("EM_Ncell",Ncell,"_Nmed",med_num_cell,"_",i_seed)
filename2 <- paste0("TCASig_Ncell",Ncell,"_Nmed",med_num_cell,"_",i_seed)
filename3 <- paste0("TCAAvgRank_Ncell",Ncell,"_Nmed",med_num_cell,"_",i_seed)

write.table(total.sig.em, filename1, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(total.sig.tca, filename2, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(avg.sig.tca, filename3, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

#write.table(sig.list, filename2, append=TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
#capture.output(summary(sig.list), file = paste0("myEMtestrun",i_seed,".txt"))
#write.table(as.data.frame(sig.list),file="myEMtest.csv", quote=F,sep=",",row.names=F)


