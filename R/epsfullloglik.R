# Log-likelihood function for EPS-full
# No interactions, no confounding

epsfullloglik_xe = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,n_cc,n_ic,
                            genotypes,ng){
    ###############################################################
    # Environmental covariates (xe) present
    ###############################################################

    len = len = 1 + dim(x_cc)[2] + dim(g_cc)[2]

    m = parameters[(len+2):length(parameters)]
    param = parameters[1:(len+1)]

    nug = length(m)+1
    probs = c()
    probs[1] = exp(m[1])/(1+exp(m[1]))
    for(k in 2:(nug-1)){
        probs[k] = (1-sum(probs))*exp(m[k])/(1+exp(m[k]))
    }

    genoprobs = c(probs,(1-sum(probs)))  # Obs: sometimes last prob less than zero, and then warning.

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp +
            exp(-(0.5*(y_ic - alpha - x_ic%*%betaE -
                           c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs

    ind = unlist(apply(g_cc,1,getgenoprob,genotypes))
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}

epsfullloglik = function(parameters,y_cc,y_ic,g_cc,n_cc,n_ic,
                         genotypes,ng){
    ###############################################################
    # Environmental covariates not present
    ###############################################################

    len = 1+dim(g_cc)[2]

    m = parameters[(len+2):length(parameters)]
    param = parameters[1:(len+1)]

    nug = length(m)+1
    probs = c()
    probs[1] = exp(m[1])/(1+exp(m[1]))
    for(k in 2:(nug-1)){
        probs[k] = (1-sum(probs))*exp(m[k])/(1+exp(m[k]))
    }

    genoprobs = c(probs,(1-sum(probs)+0.001))

    lenp = length(param)
    alpha = param[1]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp +
            exp(-(0.5*(y_ic - alpha -
                           c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    ind = unlist(apply(g_cc,1,getgenoprob,genotypes))
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}


