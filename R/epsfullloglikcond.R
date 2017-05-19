# Log-likelihood function for EPS-full
# Confounding assumed, no interactions

epsfullloglikcond = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,n_cc,n_ic,
                             genotypes,ng,xu,nu,cind){

    len = len = 1 + dim(x_cc)[2] + dim(g_cc)[2]

    m = parameters[(len+2):length(parameters)]
    param = parameters[1:(len+1)]

    allgenoprobs = list()

    # For each unique level of xe, each genetic covariate has a
    # probability distribution, so that each unique combination of
    # genotypes has a probability. These are stored in allgenotypes and
    # allgenoprobs

    for(i in 1:nu){
        # For each level of xe, get probabilities
        # for each combination of genotypes

        mtemp = m[1:(dim(genotypes)[1]-1)]
        m = m[dim(genotypes)[1]:length(m)]

        nug = length(mtemp)+1
        probs = c()
        probs[1] = exp(mtemp[1])/(1+exp(mtemp[1]))
        for(k in 2:(nug-1)){
            probs[k] = (1-sum(probs))*exp(mtemp[k])/(1+exp(mtemp[k]))
        }
        allgenoprobs[[i]] = c(probs,1-sum(probs))
    }

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    temp = 0
    temp2 = 0

    for(k in 1:nu){
        # Add up contribution to log-likelihood for all individuals
        # with missing genotypes, summing over all combinations of genotypes
        # Blocks of individuals with the same level of xe are
        # considered simultaneously
        idrows = colSums((data.frame(t(x_ic[,cind])) == t(xu[k,]))) # OBS x_ic
        id = idrows > 0 # index for individuals in block k
        genoprobs = allgenoprobs[[k]]
        tempk = 0
        for (j in 1:dim(genotypes)[1]){
            tempk = tempk +
                exp(-(0.5*(y_ic[id] - alpha - as.matrix(x_ic[id,])%*%betaE - c(t(genotypes[j,])%*%betaG))^2/sigma2))*genoprobs[j]
        }
        temp = temp + sum(log(tempk*(1/(sqrt(2*pi)*sigma))))

        # Add up contribution of completely observed individuals
        # with respect to the log of probability of their individual genotypes.
        idrows = colSums((data.frame(t(x_cc[,cind])) == t(xu[k,]))) # OBS x_cc
        id = idrows > 0 # index for individuals in block k
        ind = unlist(apply(as.matrix(g_cc[id,]),1,getgenoprob,genotypes))

        temp2 = temp2 + sum(log(genoprobs[ind]))
    }

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) - sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}
