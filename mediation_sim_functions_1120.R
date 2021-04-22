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
                              med.exp.M, med.sites, nonmed.sites, 
                              tauvec, med = TRUE){
  N_feature = dim(pure_base)[1] # number of features
  L = dim(pure_base)[2] # number of pure tissues

  tissue = matrix(0,N_feature,L)
  N_med_feature = length(med.sites)
  N_nonmed_feature = length(nonmed.sites)

  ## first generate all data, without considering mediation
  for(i in 1:L) {
      param = betaMOM(pure_base[,i],pure_sd[,i]^2)
      tissue[,i] = rbeta(N_feature,param[,1],param[,2])+rnorm(N_feature,0,tauvec)
  }

  if (med){
      ## multiply med.exp.M for mediation sites
      tissue[med.sites,] = tissue[med.sites,] * med.exp.M + rnorm(Ncell,0,tauvec)

  }

  tissue[tissue < 0] = 0.01
  tissue[tissue > 1] = 0.99

  return(tissue)
}


getProportion <- function(N_sample, cc = 100, E.exp, E.unexp, med.pi,
                          nonmed.pi, med_exp_pi, med = TRUE, Ncell){
  require("gtools")
  
  alpha.ctr <- if(Ncell==4) {c(0.1, 0.2, 0.3, 0.4)} else if(Ncell==3){
              c(0.2, 0.3, 0.5)} else {c(0.4, 0.6)}

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
getSampleMix <- function(N_sample, pure_base, pure_sd, sigma, tauvec,
                         E.exp, E.unexp, med.pi, nonmed.pi, med_exp_pi,
                         med.sites, nonmed.sites, med_exp_M,
                         is_pi_med, is_M_med){
  K = ncol(pure_base) ## number of cell types
  p = nrow(pure_base) ## number of CpG sites

  ## get proportions
  trueProp = getProportion(N_sample, cc=100, E.exp, E.unexp, med.pi,
                           nonmed.pi, med_exp_pi, med=is_pi_med, K)
  alltmp = matrix(0, p, N_sample*K)

  ## get mix
  obs.Y = matrix(0, p, N_sample)

  for(n in E.exp){ #exposed group  - relationship with mediator
    tmp = getOnePureRefPanel(pure_base, pure_sd, med.pi, nonmed.pi,
                             med.exp.M, med.sites, nonmed.sites, tauvec,
                             med = is_M_med)
    #tmp = tmp + matrix(rnorm(p*K,0,noise_sd),p,K)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, sigma)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }
  for(n in E.unexp){ #unexposed group - no relationship with mediator
    tmp = getOnePureRefPanel(pure_base, pure_sd, med.pi, nonmed.pi,
                             med.exp.M, med.sites, nonmed.sites, tauvec,
                             med = FALSE)
    #tmp = tmp + matrix(rnorm(p*K,0,noise_sd),p,K)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, sigma)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }

  obs.Y[obs.Y < 0] = 0.01
  obs.Y[obs.Y > 1] = 0.99
  rownames(obs.Y) = rownames(alltmp) = rownames(tmp) = rownames(pure_base)
  return(list(obs.Y = obs.Y, trueProp = trueProp,trueMethy = alltmp,trueMethy1 = tmp))
}



## EM Function with bootstrap - Julia
julia_command("
    function IndefBoot(M, Y, E, beta0, beta1, prop, tausq,
                      sigmasq, gammasq, theta0, theta1,
                      theta2; b=1000)

        #= ========================= =#
        #= Define loops as functions =#
        #= ========================= =#

        function lobsfun(Nsample,Ymean,Yvar,Y,Mmean,Mvar,M)
            lobsf = 0
            for i in 1:Nsample
              ypart = logpdf(Normal(Ymean[i],sqrt(Yvar)),Y[i])
              mpart = logpdf(Normal(Mmean[i],sqrt(Mvar[i])),M[i])
              lobsf += ypart + mpart
            end
            return lobsf
        end

        function calcSigmastar(Nsample,Ncell,prop,sigmasq,Sigmastarinvterm1,Sigmastarinvterm3)
            Sigmastar = Array{Float64}(undef,Ncell,Ncell,Nsample)
            for i in 1:Nsample
                  Sigmastarinvterm2 = prop[i,:]*prop[i,:]' / sigmasq
                  Sigmastar[:,:,i] = inv(Sigmastarinvterm1 + Sigmastarinvterm2 + Sigmastarinvterm3)
            end
            return Sigmastar
        end

        function calcmustarsig(Nsample,Ncell,Y,theta0,theta1,theta2,gammasq,M,prop,
                                sigmasq,beta0,beta1,E,tausq)
              mustarsig = Array{Float64}(undef,Nsample,Ncell)
              for i in 1:Nsample
                  term1 = ((Y[i] .- theta0 .- theta1 .* E[i]) .* theta2) / gammasq
                  term1 = copy(term1')
                  term2 = M[i] * prop[i,:] / sigmasq
                  term2 = copy(term2')
                  term3 = copy((beta0 + beta1 * E[i])')*inv(diagm(vec(tausq')))
                  mustarsig[i,:] = term1 + term2 + term3
              end
              return mustarsig
        end

        function calcmustar(Nsample,Ncell,mustarsig,Sigmastar)
            mustar = Array{Float64}(undef,Nsample,Ncell)
            for i in 1:Nsample
                  mustar[i,:] = mustarsig[i,:]'*Sigmastar[:,:,i]
            end
            return mustar
        end

        function calctheta2denom(Nsample,Ncell,Sigmastar,mustar)
            theta2denom = Array{Float64}(undef,Ncell,Ncell,Nsample)
            for i in 1:Nsample
                theta2denom[:,:,i] = Sigmastar[:,:,i] .+ mustar[i,:]*mustar[i,:]'
            end
            return theta2denom
        end

        function calcsigmavec(Nsample,prop,Sigmastar,M,mustar)
            sigmavec = Array{Float64}(undef,Nsample)
            for i in 1:Nsample
                sigmavec[i] = prop[i,:]'*Sigmastar[:,:,i]*prop[i,:] .+ (M[i] .-
                mustar[i,:]'*prop[i,:]) .^2
            end
            return sigmavec
        end

        function calcgammavec(theta2new,Sigmastar,Y,theta0new,theta1new,mustar,E)
            gammavec = Array{Float64}(undef,Nsample)
            for i in 1:Nsample
                gammavec[i] = theta2new'*Sigmastar[:,:,i]*theta2new .+
                        (Y[i] - theta0new - theta1new*E[i] - (mustar[i,:]'*theta2new)[1]) ^2
            end
            return gammavec
        end

        function calctauvec(Nsample,Ncell,E,mustar,beta0new,beta1new,Sigmastar)
            tauvec = Array{Float64}(undef,Nsample,Ncell)
            for i in 1:Nsample, j in 1:Ncell
                tauvec[i,j] = (mustar[i,j] - beta0new[j] - beta1new[j]*E[i])^2+
                                      Sigmastar[:,:,i][j,j]
            end
            return tauvec
        end



        #= ========================= =#
        #= Define EM Function        =#
        #= ========================= =#

                function CTMEM(M, Y, E, beta0, beta1, prop, tausq,
                              sigmasq, gammasq, theta0, theta1,
                              theta2; maxiter=500, tol=0.001)
                iter = 1
                err = 2000
                Nsample = length(Y)
                Ncell = size(prop)[2]
                Emat = repeat(E,1,Ncell)

                          #= ========================== =#
                          #= Define likelihood Function =#
                          #= ========================== =#

                          function calc_lobs(M, Y, E, beta0, beta1, prop, tausq,
                                sigmasq, gammasq, theta0, theta1,theta2)

                                Nsample = length(Y)
                                Ncell = size(prop)[2]
                                Emat = repeat(E,1,Ncell)

                                # calculate observed likelihood terms
                                beta1E = Emat .* beta1'
                                beta0beta1E = beta1E .+ beta0'
                                theta2beta = beta0beta1E .* theta2'

                                Ymean = theta0 .+ theta1 .* E .+ sum(theta2beta,dims=2)
                                Yvar = sum(theta2' .^2 .* tausq') + gammasq

                                Mmean = sum(beta0beta1E .* prop,dims=2)
                                Mvar = sum(prop .^ 2 .* tausq',dims=2) .+ sigmasq


                                lobs = lobsfun(Nsample,Ymean,Yvar,Y,Mmean,Mvar,M)

                                return lobs

                          end  ### end likelihood function

                # initialize likelihood for EM
                lobsi = calc_lobs(M, Y, E, beta0, beta1, prop, tausq,
                       sigmasq, gammasq, theta0, theta1, theta2)


                while iter<=maxiter && err>=tol

                    #= ======= =#
                    #= E step  =#
                    #= ======= =#

                    #calculate mu* and sigma* (mean and variance of conditional distribution of m_k)
                      Sigmastarinvterm1 = theta2*theta2' / gammasq
                      Sigmastarinvterm3 = inv(diagm(vec(tausq')))
                      Sigmastar = calcSigmastar(Nsample,Ncell,prop,sigmasq,Sigmastarinvterm1,
                                                                      Sigmastarinvterm3)

                      mustarsig = calcmustarsig(Nsample,Ncell,Y,theta0,theta1,theta2,gammasq,M,prop,
                                sigmasq,beta0,beta1,E,tausq)

                      mustar = calcmustar(Nsample,Ncell,mustarsig,Sigmastar)

                    #= ======= =#
                    #= M step  =#
                    #= ======= =#

                    # beta0 update
                    beta1E = Emat .* beta1'
                    beta0new = sum(mustar-beta1E,dims=1) / Nsample
                    beta0new = dropdims(beta0new';dims=2)


                    # beta1 update
                    beta0mat = repeat(beta0new',Nsample,1)
                    beta1new = sum((mustar-beta0mat) .* Emat, dims=1) ./ sum(Emat .^ 2,dims=1)
                    beta1new = dropdims(beta1new';dims=2)

                    # theta0 update
                    mutheta2 = mustar .* theta2'
                    theta0new = sum(Y .- theta1 .* E - sum(mutheta2,dims=2)) / Nsample

                    # theta1 update
                    theta1new = sum((Y .- theta0new - sum(mutheta2,dims=2)) .* E) / sum(E .^ 2)

                    # theta2 update
                      # denominator
                      theta2denom = calctheta2denom(Nsample,Ncell,Sigmastar,mustar)
                      denom = inv(dropdims(sum(theta2denom,dims=3);dims=3))
                    numer = (Y .- theta0new .- theta1new .* E) .* mustar
                    theta2new = sum(numer,dims=1)*denom
                    theta2new = dropdims(theta2new';dims=2)


                    # sigma^2 update
                    sigmavec = calcsigmavec(Nsample,prop,Sigmastar,M,mustar)
                    sigmasqnew = sum(sigmavec) / Nsample

                    # gamma^2 update
                    gammavec = calcgammavec(theta2new,Sigmastar,Y,theta0new,theta1new,mustar,E)
                    gammasqnew = sum(gammavec) / Nsample

                    # tau^2 update
                    tauvec = calctauvec(Nsample,Ncell,E,mustar,beta0new,beta1new,Sigmastar)
                    tausqnew = sum(tauvec,dims=1) / Nsample
                    tausqnew = dropdims(tausqnew';dims=2)


                    # Evaluate new observed log-likelihood
                        lobsnew = calc_lobs(M, Y, E, beta0new, beta1new, prop,
                                       tausqnew, sigmasqnew, gammasqnew,
                                       theta0new, theta1new, theta2new)

                err = abs(lobsnew - lobsi)

                #Core.println(lobsnew)

                #set new params
                beta0 = beta0new
                beta1 = beta1new
                theta0 = theta0new
                theta1 = theta1new
                theta2 = theta2new
                tausq = tausqnew
                sigmasq = sigmasqnew
                gammasq = gammasqnew

                lobsi = lobsnew

                iter += 1

                end



                return beta0, beta1, theta0, theta1, theta2, tausq, sigmasq, gammasq, iter

                end # End EM Function

        #= ========================= =#
        #= Define Bootstrap Function =#
        #= ========================= =#

        function calcboot(b,E,Y,M,Nsample,Ncell, beta0, beta1, prop, tausq,
                sigmasq, gammasq, theta0, theta1, theta2)
            indefb::SharedArray{Float64,2}=zeros(b,Ncell)
            @sync @distributed for i in 1:b

              theIndex = sample(1:Nsample,Nsample,replace=true)

              Eb = E[theIndex]
              Yb = Y[theIndex]
              Mb = M[theIndex]
              propb = prop[theIndex,:]

              EMboot = CTMEM(Mb, Yb, Eb, beta0, beta1, propb, tausq,
                sigmasq, gammasq, theta0, theta1, theta2)

              indefb[i,:] = EMboot[2] .* EMboot[5]
            end
            
            #display(indefb)
            return indefb
        end



        Ncell = size(prop)[2]
        Nsample = size(prop)[1]

        myindefb = calcboot(b,E,Y,M,Nsample,Ncell, beta0, beta1, prop, tausq,
                sigmasq, gammasq, theta0, theta1, theta2)

        EMout = CTMEM(M, Y, E, beta0, beta1, prop, tausq,
               sigmasq, gammasq, theta0, theta1, theta2)

        ObsIE = EMout[2] .* EMout[5]

        function getquant(indefb,Ncell)
            bootint = Array{Float64}(undef,2,Ncell)
            pvals = Array{Float64}(undef,1,Ncell)
            for i in 1:Ncell
                bootint[:,i] = quantile(myindefb[:,i],[.025/Ncell,1-(.025/Ncell)])
                pvals[i] = 2*(min(length(myindefb[:,i][myindefb[:,i] .> 0]), 
                length(myindefb[:,i][myindefb[:,i] .< 0]))/length(myindefb[:,i]))
            end
            return bootint, pvals
        end

        mybootint = getquant(myindefb,Ncell)[1]
        mybootpvals = getquant(myindefb,Ncell)[2]

        function calcindefsig(Ncell,bootint)
            indefSig = Array{Float64}(undef,Ncell)
            for i in 1:Ncell
                if mybootint[1,i] < 0 && mybootint[2,i] > 0
                  indefSig[i] = 0
                else
                  indefSig[i] = 1
                end
            end
            return indefSig
        end

        indefSig = calcindefsig(Ncell,mybootint)

        return indefSig, ObsIE, mybootint, mybootpvals


end")

julia_command("
    function EMBoot(M, Y, E, coefs, prop, medpi, initmethod, thetaEffect, toasttheta2)


    Ncpg = size(coefs)[2]
    Ncell = size(prop)[2]
    BootSig = Array{Float64}(undef,Ncpg,Ncell)
    Bootpvals = Array{Float64}(undef,Ncpg,Ncell)

    # set initial values

    tausq = repeat(rand(InverseGamma(5,0.004),1),Ncell)
    sigmasq = rand(InverseGamma(5,0.004),1)[1]
    gammasq = rand(InverseGamma(5,0.004),1)[1]
    theta0 = 0
    theta1 = 0
    
    if initmethod == \"rand\"
        theta2 = rand(Normal(0.3,0.05),Ncell)
      elseif initmethod == \"true\"
        theta2 = zeros(Ncell)
        medpi = Int8(medpi)
        theta2[medpi] = thetaEffect
      else
        theta2 = toasttheta2
    end

    # loop through CpG sites
    for i in 1:Ncpg
          Core.println(i)

          beta0 = coefs[1:Ncell,i]
          beta1 = coefs[(Ncell+1):(2*Ncell),i]
           
          myres = IndefBoot(M[i,:], Y, E, beta0, beta1, prop, tausq,
                       sigmasq, gammasq, theta0, theta1, theta2)
          BootSig[i,:] = myres[1]
          Bootpvals[i,:] = myres[4]
    end

    return BootSig, Bootpvals

end")


getIE <- function(O,E,est.Mik){
  Ncell <- length(est.Mik)
  med_num_M <- dim(est.Mik[[1]])[1]
  Nsample <- dim(est.Mik[[1]])[2]
  
  Mikmat <- matrix(unlist(est.Mik),nrow=Nsample,ncol=Ncell)
  
  obsIE <- matrix(NA,nrow=med_num_M,ncol=Ncell) #store observed indirect effect for each site
  for(i in 1:med_num_M){
    betas <- matrix(NA,Ncell,1)
    for(j in 1:Ncell){
      medmod <- summary(fastLm(est.Mik[[j]][i,]~E))
      betas[j,] <- medmod$coefficients[2,1]
    }
    
    outmod <- summary(fastLm(O~E+Mikmat[,seq(i,Ncell*med_num_M,by=med_num_M)]))
    theta2 <- matrix(outmod$coefficients[3:(Ncell+2),1],nrow=1,ncol=Ncell)
    obsIE[i,] <- betas*t(theta2)
  }
  return(list(betas=betas,theta2=theta2,obsIE=obsIE))
}

getTOAST <- function(E,O,M,prop){
  design <- data.frame(E = as.factor(E))
  Design_out <- makeDesign(design, prop)
  fm <- fitModel(Design_out, M)
  allbetas <- fm$coefs
  betas <- fm$coefs[(Ncell+1):(2*Ncell),] ## dim: 2*Ncell X NCpG
  
  design2 <- data.frame(O=O,E = as.factor(E))
  Design_out2 <- makeDesign(design2, prop)
  fm2 <- fitModel(Design_out2, M)
  thetas <- fm2$coefs[(Ncell+1):(2*Ncell),] ## dim: 2*Ncell X NCpG
  
  indef <- betas*thetas
  
  return(list(allbetas=allbetas,betas=betas,thetas=thetas,indef=indef))
}

