library(TOAST)
load("testTOAST.rda")
source("fitModel.R")

design <- data.frame(E = as.factor(E))
Prop <- trueProp
colnames(Prop) <- c("CD8T","CD4T","NK","Bcell") ## made up because need colnames
Design_out <- makeDesign(design, Prop)

Design_out$design_matrix ## why is there only 7 columns? Do we have intercept in the model? I thought we don't have intercept?


fm <- fitModel(Design_out, Y.raw)

coefs <- fm$coefs ## dim: 2*Ncell X NCpG

fm$all_coefs
fm$all_cell_types

### following two lines don't work
res <- csTest(fm, coef="E", sort = FALSE)
res <- csTest(fm, coef="E", cell_type="CD8T", sort = FALSE)


##########################
load("testTOAST.rda")

## check the true cell type specific methylation
## med.sites should have cs methylation correlated with E
par(mfrow=c(2,2), ask=TRUE)
for( i in 1:nrow(trueMethy) ) {
    trueMethy.matrix = matrix(trueMethy[i,], ncol=4)
    for(j in 1:ncol(trueMethy.matrix))
        boxplot(trueMethy.matrix[,j] ~ E)
}

## Create mixing data
p = nrow(trueMethy)
N_sample = 100
obs.Y = matrix(0, p, N_sample)
noise_sd = 0.01
for(n in 1:p) {
    trueMethy.ind = matrix(trueMethy[n,], ncol=4)
    obs.Y[n,] = trueMethy.ind %*% trueProp[n,] + rnorm(N_sample, 0, noise_sd)
}
obs.Y[obs.Y < 0] = 0.01
obs.Y[obs.Y > 1] = 0.99
rownames(obs.Y) = rownames(trueMethy)

#### Doing TOAST test
Y.raw = obs.Y
design <- data.frame(E = as.factor(E))
Prop <- trueProp
colnames(Prop) <- c("CD8T","CD4T","NK","Bcell") ## made up because need colnames
Design_out <- makeDesign(design, Prop)

fm <- fitModel(Design_out, Y.raw)

coefs <- fm$coefs

fm$all_coefs
fm$all_cell_types

### following two lines don't work
res <- csTest(fm, coef="E", sort = FALSE)

## med.sites is 4 and 8. The 8th cg is not significant, why?

boxplot(trueProp)

## Maybe because the proportion is too low, and sample size is small.
