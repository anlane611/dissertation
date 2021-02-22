source("sim_functions.R")

load("BloodPureDNAmethy_4ct.rda")
pure.sd[pure.sd == 0] = 1e-6 ## avoid NA for estimating parameters of beta distribution
Ncell = 4 ## number of cell types
Nsample = 50 ## number of samples

out = getSampleMix(N_sample = Nsample, 
                   pure.mean[,1:Ncell], pure.sd[,1:Ncell], 
                   noise_sd = 0.01)
Y.raw = out$obs.Y ## observed DNA methylation samples
trueProp = out$trueProp ## proportions in samples