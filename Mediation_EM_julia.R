#### EM Algorithm for Cell-type-specific mediation

library(JuliaCall)
julia <- JuliaCall::julia_setup("/Applications/JuliaPro-1.5.2-1.app/Contents/Resources/julia/Contents/Resources/julia/bin")
julia_library("LinearAlgebra")
julia_library("Distributions")

library(invgamma)

## initial values pulled from simulation code

## mediating sites 5, 9; mediating cell types 1 and 4 respectively
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


#testing initial values
theta0 <- 0
theta1 <- 0


tausq <- rep(rinvgamma(1,5,0.04),Ncell)
sigmasq <- rinvgamma(1,5,0.04)
gammasq <- (sd(O3)/5)^2



theta2 <- c(0,0,0,0) #case 1
site9cell1.lowvar.1 <- julia_call("IndefBoot",M, Y, E, beta0, beta1, pi, tausq,
                                  sigmasq, gammasq, theta0, theta1, theta2)

theta2 <- c(0.2,0.2,0.2,0.2) #case 2
site9cell1.lowvar.2 <- julia_call("IndefBoot",M, Y, E, beta0, beta1, pi, tausq,
                                  sigmasq, gammasq, theta0, theta1, theta2)

theta2 <- c(0.3,0.1,0.1,0.1) #case 3
site9cell1.lowvar.3 <- julia_call("IndefBoot",M, Y, E, beta0, beta1, pi, tausq,
                                  sigmasq, gammasq, theta0, theta1, theta2)

theta2 <- c(0.4,0,0,0) #case 4
site9cell1.lowvar.4 <- julia_call("IndefBoot",M, Y, E, beta0, beta1, pi, tausq,
                                  sigmasq, gammasq, theta0, theta1, theta2)

site9cell1.lowvar <- rbind(site9cell1.lowvar.1,site9cell1.lowvar.2,site9cell1.lowvar.3,site9cell1.lowvar.4)


#true values
theta0 <- 0
theta1 <- 0

theta2 <- c(0,0,0,0.4)

tausq <- rep((0.1)^2,Ncell)
sigmasq <- 0.1^2
gammasq <- (sd(O3)/5)^2



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
                lobsvec = Array{Float64}(undef,maxiter,6)
                
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
                lobsvec[iter,6] = iter
                
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
      
        mylobsvec = EMout[5]

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

