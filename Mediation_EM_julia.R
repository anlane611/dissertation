#### EM Algorithm for Cell-type-specific mediation

library(JuliaCall)
julia <- JuliaCall::julia_setup("/Applications/JuliaPro-1.5.2-1.app/Contents/Resources/julia/Contents/Resources/julia/bin")
julia_library("LinearAlgebra")
julia_library("Distributions")

library(invgamma)

## initial values pulled from simulation code
load("Nsample10000.RData")
load("Nsample500.Rdata")

M <- Y.raw
Y <- O
E <- E

#get one cpg site
i <- 9
M <- M[i,]
beta0 <- coefs[1:Ncell,i]
beta1 <- coefs[(Ncell+1):(2*Ncell),i]
pi <- trueProp
# theta2 <- theta.2[i,]


# theta0 <- theta0[i]
# theta1 <- theta1[i]
truem <- matrix(trueMethy[i,],Nsample,Ncell)

#true values
theta0 <- 0
theta1 <- 0
theta2 <- c(0,0,0,0.4)

tausq <- rep((0.1)^2,Ncell)
sigmasq <- 0.1^2
gammasq <- (sd(O3)/5)^2

## use true values as initial values
res <- CTMEM(M=M, Y=Y, E=E, beta0=beta0, beta1=beta1, pi=pi, tausq=tausq,
             sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
             theta2=theta2, max.iter=200, tol=0.001)



## function to calculate observed likelihood
calc.lobs <- function(M, Y, E, beta0, beta1, pi, tausq,
                    sigmasq, gammasq, theta0, theta1,
                    theta2){

  Nsample <- length(Y)
  Ncell <- ncol(pi)
  E.mat <- matrix(E,nrow=Nsample,ncol=Ncell)

  ## initialize observed log-likelihood
  beta1E <- sweep(E.mat,2,beta1,FUN="*") #beta1^kE
  beta0beta1E <- sweep(beta1E,2,beta0,FUN="+") #beta0^k + beta1^kE
  theta2beta <- sweep(beta0beta1E,2,theta2,FUN="*") #theta2^k + (beta0^k + beta1^kE)
  
  Y.mean <- theta0 + theta1*E + rowSums(theta2beta) #sum_k to get mean of Y
  Y.var <- sum(theta2^2*tausq)+gammasq #sum_k theta2^2(k)*tau_k^2 + gamma^2
  
  M.mean <- rowSums(beta0beta1E*pi) #sum_k pi_k*(beta0^k+beta1^kE)
  M.var <- rowSums(sweep(pi^2,2,tausq,FUN="*")) + sigmasq #sum_k pi_k^2*tau_k^2 + sigma^2
  
  l.obs <- 0
  for(i in 1:Nsample){
    y.part <- log(dnorm(Y[i],Y.mean[i],sqrt(Y.var)))
    m.part <- log(dnorm(M[i],M.mean[i],sqrt(M.var[i])))
    l.obs <- l.obs + y.part + m.part
  }

  return(l.obs)

}


#input initial values (from TCA) for now
CTMEM <- function(M, Y, E, beta0, beta1, pi, tausq,
                  sigmasq, gammasq, theta0, theta1,
                  theta2, max.iter=1000, tol=0.001){
  err <- 1000
  iter <- 1
  Nsample <- length(Y)
  Ncell <- ncol(pi)
  E.mat <- matrix(E,nrow=Nsample,ncol=Ncell)
  
  ## initialize observed log-likelihood
  l.obs <- calc.lobs(M=M, Y=Y, E=E, beta0=beta0, beta1=beta1, pi=pi, tausq=tausq,
                     sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
                     theta2=theta2)
  
  ## initialize some temporary variables. Do this outside the EM loop
  sigma.vec <- gamma.vec <- numeric(Nsample)
  Sigmastar <- matrix(0, Ncell, Ncell*Nsample)
  mustar.sig <- matrix(0,Nsample,Ncell)
  mustar <- matrix(0,Nsample,Ncell)
  theta2.denom <- vector("list", Nsample)
  tau.vec <- matrix(0,Nsample,Ncell)
  
  while(err>=tol & iter<=max.iter){
    ## E step: calculate mu* and sigma* (mean and variance of conditional distribution of m_k)
    
    Sigmastarinv.term1 <- (theta2%*%t(theta2))/gammasq
    Sigmastarinv.term3 <- solve(diag(tausq))
    for(i in 1:Nsample){
      Sigmastarinv.term2 <- (pi[i,]%*%t(pi[i,]))/sigmasq
      Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)] <- solve(Sigmastarinv.term1+
                                                         Sigmastarinv.term2+Sigmastarinv.term3)
    }
    
    for(i in 1:Nsample){
      term1 <- ((Y[i]-theta0-theta1*E[i])*t(theta2))/gammasq
      term2 <- (M[i]*t(pi[i,]))/sigmasq
      term3 <- t(beta0+beta1*E[i])%*%solve(diag(tausq))
      mustar.sig[i,] <- term1 + term2 + term3
    }
    
    for(i in 1:Nsample){
      mustar[i,] <- mustar.sig[i,]%*%Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)]
    }
    
    
    ## M-step
    
    # update beta0
    beta1E <- sweep(E.mat,2,beta1,FUN = "*")
    beta0.new <- (colSums(mustar-beta1E))/Nsample
    
    # update beta1
    beta0mat <- matrix(beta0.new,nrow=Nsample,ncol=Ncell,byrow = TRUE)
    beta1.new <- (colSums((mustar-beta0mat)*E.mat))/colSums(E.mat^2)
    
    # update theta0
    mu.theta2 <- sweep(mustar,2,theta2,FUN="*")
    theta0.new <- sum(Y-theta1*E-rowSums(mu.theta2))/Nsample
    
    # update theta1
    theta1.new <- sum((Y-theta0.new-rowSums(mu.theta2))*E)/sum(E^2)
    
    # update theta2 (k-vector)
    for(i in 1:Nsample){ #calculate denominator for each i
      theta2.denom[[i]] <- Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)]+
        (mustar[i,]%*%t(mustar[i,]))
    }
    denom <- solve(Reduce("+",theta2.denom)) #sum over i and take inverse
    numer <- sweep(mustar,1,Y-theta0.new-theta1.new*E,FUN="*") #numerator: (Y-theta0-theta1E)mu*
    theta2.new <- colSums(numer)%*%denom #sum numerator over i and multiply by denominator
    theta2.new <- t(theta2.new)
    
    # update sigma^2
    for(i in 1:Nsample){
      sigma.vec[i] <- t(pi[i,])%*%Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)]%*%pi[i,] +
        (M[i] - t(mustar[i,])%*%pi[i,])^2
    }
    sigmasq.new <- sum(sigma.vec)/Nsample
    
    # update gamma^2
    for(i in 1:Nsample){
      gamma.vec[i] <- t(theta2.new)%*%Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)]%*%theta2.new+
        (Y[i] - theta0.new - theta1.new*E[i] - t(mustar[i,])%*%theta2.new)^2
    }
    gammasq.new <- sum(gamma.vec)/Nsample
    
    # update tau^2 (k-vector)
    for(i in 1:Nsample){
      for(j in 1:Ncell){
        tau.vec[i,j] <- (mustar[i,j] - beta0.new[j] - beta1.new[j]*E[i])^2+
          Sigmastar[,(Ncell*i-Ncell+1):(Ncell*i)][j,j]
      }
    }
    tausq.new <- colSums(tau.vec)/Nsample
    
    ## evaluate observed log-likelihood to determine convergence
    l.obs.new <- calc.lobs(M=M, Y=Y, E=E, beta0=beta0.new, beta1=beta1.new, pi=pi,
                           tausq=tausq.new, sigmasq=sigmasq.new, gammasq=gammasq.new,
                           theta0=theta0.new, theta1=theta1.new, theta2=theta2.new)
    
    err <- abs(l.obs.new - l.obs)
    
    # err2 <- sum(beta0.new-beta0)^2 + sum(beta1.new-beta0)^2 + (theta0.new-theta0)^2 +
    #         (theta1.new-theta1)^2 + sum(theta2.new-theta2)^2 + (gammasq.new-gammasq)^2 +
    #         (sigmasq.new-sigmasq)^2 + sum(tausq.new-tausq)^2
    
    #cat("iter", iter, ": llik=", l.obs.new, ", err=", err, "\n")
    
    #set new params
    beta0 <- beta0.new
    beta1 <- beta1.new
    theta0 <- theta0.new
    theta1 <- theta1.new
    theta2 <- theta2.new
    tausq <- tausq.new
    sigmasq <- sigmasq.new
    gammasq <- gammasq.new
    
    l.obs <- l.obs.new
    
    iter <- iter+1
  }
  
  return(list(sigmasq=sigmasq, gammasq=gammasq, tausq=tausq,
              beta0=beta0, beta1=beta1,
              theta0=theta0, theta1=theta1, theta2=theta2,
              err=err, mustar=mustar, iter=iter))
  
}

## use true values as initial values
time0 <- Sys.time()
res <- CTMEM(M=M, Y=Y, E=E, beta0=beta0, beta1=beta1, pi=pi, tausq=tausq,
      sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
      theta2=theta2, max.iter=200, tol=0.01)
diff2 <- Sys.time()-time0


## use some random initial values
Ncell <- ncol(pi)
beta0.0 = rep(0.2, Ncell)
beta1.0 = rep(0, Ncell)
res = CTMEM(M=M, Y=Y, E=E, beta0=beta0.0, beta1=beta1.0, pi=pi, tausq=tausq,
        sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
        theta2=theta2, max.iter=1000, tol=0.01)


### bootstrap for indirect effect

IndefBoot <- function(M, Y, E, beta0, beta1, pi, tausq,
                        sigmasq, gammasq, theta0, theta1,
                        theta2, theIndex, b=100){
    Ncell <- dim(pi)[2]
    Nsample <- dim(pi)[1]
    indef.b <- matrix(0,b,Ncell) # indirect effects for bootstrap
    Indef.Sig <- c()
    

    EMout <- CTMEM(M=M, Y=Y, E=E, beta0=beta0, beta1=beta1, pi=pi, tausq=tausq,
                   sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
                   theta2=theta2, max.iter=500, tol=0.001)
    obs.indef <- EMout$beta1*EMout$theta2
    
    for(nboot in 1:b){
        print(nboot)
        # get indices for bootstrap sample
        #index.b <- sample(1:Nsample,Nsample,replace = TRUE)
        index.b <- theIndex
        E.b <- E[index.b]
        M.b <- M[index.b]
        pi.b <- pi[index.b,]
        Y.b <- Y[index.b]

        EMout.b <- CTMEM(M=M.b, Y=Y.b, E=E.b, beta0=beta0, beta1=beta1, pi=pi.b, tausq=tausq,
                 sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
                 theta2=theta2, max.iter=500, tol=0.001)

  
        indef.b[nboot,] <- EMout.b$beta1*EMout.b$theta2

    } # end b

    bootint <- apply(indef.b,2,function(x) quantile(x, c(.05/4,1-(.05/4))))

    ## indirect effect in bootstrap interval?
    for(i in 1:Ncell){
        Indef.Sig[i] <- ifelse(0 < bootint[2,i] &
                          0 > bootint[1,i],0,1)
    }
        
        return(list(Indef.Sig=Indef.Sig, indirect.effect=obs.indef, bootint=bootint, 
                    indef.b=indef.b))

}

time0 <- Sys.time()
IndefBoot(M=M, Y=Y, E=E, beta0=beta0, beta1=beta1, pi=pi, tausq=tausq,
                      sigmasq=sigmasq, gammasq=gammasq, theta0=theta0, theta1=theta1,
                      theta2=theta2, theIndex=theIndex, b=2)
diff2 <- Sys.time() - time0


### Julia ###


julia_command("
    function IndefBoot(M, Y, E, beta0, beta1, prop, tausq,
                      sigmasq, gammasq, theta0, theta1,
                      theta2; b=200)
                      
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
                              theta2; maxiter=1000, tol=0.001)
                iter = 1    
                err = 2000
                Nsample = length(Y)
                Ncell = size(prop)[2]
                Emat = repeat(E,1,4)
                lobsvec = Array{Float64}(undef,maxiter,5)
                
                          #= ========================== =#
                          #= Define likelihood Function =#
                          #= ========================== =#
                
                          function calc_lobs(M, Y, E, beta0, beta1, prop, tausq,
                                sigmasq, gammasq, theta0, theta1,theta2)
                                
                                Nsample = length(Y)
                                Ncell = size(prop)[2]
                                Emat = repeat(E,1,4)
            
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
                
                lobsvec[iter,1] = lobsnew
                lobsvec[iter,2:5] = theta2new
                
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
            
            
                    
                return beta0, beta1, theta0, theta1, theta2, tausq, sigmasq, gammasq, iter, lobsvec
                
                end # End EM Function
        
        #= ========================= =#
        #= Define Bootstrap Function =#
        #= ========================= =#

        # function calcboot(b,E,Y,M,Nsample,Ncell, beta0, beta1, prop, tausq,
        #         sigmasq, gammasq, theta0, theta1, theta2)
        #     indefb = Array{Float64}(undef,b,Ncell)
        #     for i in 1:b
        # 
        #       theIndex = sample(1:Nsample,Nsample,replace=true)
        # 
        #       Eb = E[theIndex]
        #       Yb = Y[theIndex]
        #       Mb = M[theIndex]
        #       propb = prop[theIndex,:]
        # 
        #       EMboot = CTMEM(Mb, Yb, Eb, beta0, beta1, propb, tausq,
        #         sigmasq, gammasq, theta0, theta1, theta2)
        # 
        #       indefb[i,:] = EMboot[2] .* EMboot[5]
        #     end
        #     return indefb
        # end



         Ncell = size(prop)[2]
         Nsample = size(prop)[1]

        # myindefb = calcboot(b,E,Y,M,Nsample,Ncell, beta0, beta1, prop, tausq,
        #         sigmasq, gammasq, theta0, theta1, theta2)

         EMout = CTMEM(M, Y, E, beta0, beta1, prop, tausq,
                sigmasq, gammasq, theta0, theta1, theta2)

        # ObsIE = EMout[2] .* EMout[5]
        # 
        # function getquant(indefb,Ncell)
        #     bootint = Array{Float64}(undef,2,Ncell)
        #     for i in 1:Ncell
        #         bootint[:,i] = quantile(myindefb[:,i],[.05/Ncell,1-(.05/4)])
        #     end
        #     return bootint
        # end
        # 
        # mybootint = getquant(myindefb,Ncell)
        # 
        # function calcindefsig(Ncell,bootint)
        #     indefSig = Array{Float64}(undef,Ncell)
        #     for i in 1:Ncell
        #         if mybootint[1,i] < 0 && mybootint[2,i] > 0
        #           indefSig[i] = 0
        #         else
        #           indefSig[i] = 1
        #         end
        #     end
        #     return indefSig
        # end
        # 
        # indefSig = calcindefsig(Ncell,mybootint)
        # 
        # return indefSig, ObsIE, mybootint
      
        mylobsvec = EMout[10]

        return mylobsvec
end")

julia_command("
    function EMBoot(M, Y, E, coefs, prop)
    
    # set initial values

    tausq = repeat(rand(InverseGamma(5,0.004),1),4)
    sigmasq = rand(InverseGamma(5,0.004),1)[1]
    gammasq = 0.0004351324 #(sd(O3)/5)^2
    theta0 = 0
    theta1 = 0
    theta2 = zeros(4)

    # loop through CpG sites

    Ncpg = size(coefs)[2]
    Ncell = size(prop)[2]
    BootSig = Array{Float64}(undef,Ncpg,Ncell)
    
    for i in 1:Ncpg
          Core.println(i)

          beta0 = coefs[1:Ncell,i]
          beta1 = coefs[(Ncell+1):(2*Ncell),i]
          BootSig[i,:] = IndefBoot(M[i,:], Y, E, beta0, beta1, prop, tausq,
                       sigmasq, gammasq, theta0, theta1, theta2)[1]
    end
    
    return BootSig
    
        
end")


bootint <- apply(indef.b,2,function(x) quantile(x, c(.05/4,1-(.05/4))))

time0 <- Sys.time()
test0.45 <- julia_call("IndefBoot",M, Y, E, beta0, beta1, pi, tausq,
           sigmasq, gammasq, theta0, theta1, theta2)
diff <- Sys.time()-time0


time0 <- Sys.time()
test <- julia_call("EMBoot",Y.raw[c(5,9),], O, E, coefs[,c(5,9)], Prop)
diff <- Sys.time()-time0

time0 <- Sys.time()
test <- julia_call("EMBoot",M, Y, E, coefs, pi)
diff <- Sys.time()-time0


theta1.new <- sum((Y-theta0.new-rowSums(mu.theta2))*E)/sum(E^2)

#testing initial values
beta0 <- c(0,0,0,0)
beta1 <- c(0,0,0,0)
theta0 <- 0
theta1 <- 0
theta2 <- c(-0.1,-0.1,-0.1,0.45)

tausq <- rep(rinvgamma(1,5,0.04),Ncell)
sigmasq <- rinvgamma(1,5,0.04)
gammasq <- (sd(O3)/5)^2

truevals <- c(0,0,0,0.4)
par(mfrow=c(2,2))
for(i in 2:4){
  plot(test[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("true values as initial values; max.iter=500, tol=0.01", side = 3, line = -1, outer = TRUE)


par(mfrow=c(2,2))
for(i in 2:4){
  plot(test3[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test3[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("0 as initial values; max.iter=500, tol=0.01", side = 3, line = -1, outer = TRUE)


truevals <- c(0,0,0,0.4)
par(mfrow=c(2,2))
for(i in 2:4){
  plot(test2[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test2[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("true values as initial values; max.iter=500, tol=0.0001", side = 3, line = -1, outer = TRUE)


par(mfrow=c(2,2))
for(i in 2:4){
  plot(test4[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test4[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("0 as initial values; max.iter=500, tol=0.0001", side = 3, line = -1, outer = TRUE)


truevals <- c(0,0,0,0.4)
par(mfrow=c(2,2))
for(i in 2:4){
  plot(test5[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test5[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("true values as initial values; max.iter=500, tol=0.0001", side = 3, line = -1, outer = TRUE)


par(mfrow=c(2,2))
for(i in 2:4){
  plot(test6[,i],main=paste0("Cell ",i-1),xlab="Iteration",ylab="Estimate",ylim=c(-0.4,0.3))
  abline(h=truevals[i-1],col="red")
}
plot(test6[,5],main="Cell 4 (mediator)",xlab="Iteration",ylab="Estimate",ylim=c(0.1,0.5))
abline(h=0.4,col="red")
mtext("0 as initial values; max.iter=2000, tol=0.0001", side = 3, line = -1, outer = TRUE)
