# Log-likelihood function for EPS-complete
# Interactions, no confounding

epscomploglikint = function(parameters,data,ng,interactind,
                            geneffect = "additive"){
    len = dim(data)[2]
    neg = length(interactind)
    m = parameters[(len+1+neg+1):length(parameters)]
    param = parameters[1:(len+1+neg)]

    # Extract all combinations of genotypes, and the
    # corresponding probability for each combination.
    # It is assumed that all genetic variables are statistically
    # independent. Use function genprob().
    if(geneffect == "additive"){
        probs = c()
        count = 1
        for(t in 1:ng){
            probs[count] = exp(m[count])/(1+exp(m[count]))
            probs[count + 1] = (1-probs[count])*
                exp(m[count+1])/(1+exp(m[count+1]))
            count = count + 2
        }
        getgeno = genprob(ng,probs)
        genotypes = as.matrix(getgeno[[1]])
        genoprobs = getgeno[[2]]
        #print(genoprobs)
    }else{
        probs = 0.9*exp(m)/(1+exp(m))
        getgeno = genprob(ng,probs, geneffect = geneffect)
        genotypes = as.matrix(getgeno[[1]])
        genoprobs = getgeno[[2]]
    }



    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]

    y_cc = data_cc[,1]
    y_ic = data_ic[,1]
    g_cc = as.matrix(data_cc[,(len-ng+1):len])
    x_cc = as.matrix(data_cc[,2:(len-ng)])
    x_ic = as.matrix(data_ic[,2:(len-ng)])
    n_cc = length(y_cc)
    n_ic = length(y_ic)

    neg = length(interactind)
    xg_cc = matrix(NA,ncol = length(interactind),nrow = n_cc)
    for(i in 1:neg){
        xg_cc[,i] = x_cc[,interactind[[i]][2]]*g_cc[,interactind[[i]][1]]
    }
    interactgenotypes = list()

    for (i in 1:dim(genotypes)[1]){
        interactgenotypes[[i]] = matrix(NA,ncol = length(interactind),
                                        nrow = n_ic)
        for(j in 1:neg){
            interactgenotypes[[i]][,j] =
                x_ic[,interactind[[j]][2]]*genotypes[i,interactind[[j]][1]]
        }
    }

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    betaEG = param[(len+1):(len+length(interactind))]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp +
            exp(-(0.5*(y_ic - alpha - x_ic%*%betaE -
                    c(t(genotypes[i,])%*%betaG) -
                    interactgenotypes[[i]]%*%betaEG)^2)/sigma2)*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs
    if(geneffect=="additive"){rowind = 3}else{rowind = 2}
    ind = rep(1,n_cc)
    for(j in 1:ng){
        ind = ind + g_cc[,j]*(rowind)^(j-1)
    }
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - x_cc%*%betaE -
                     g_cc%*%betaG - xg_cc%*%betaEG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}
