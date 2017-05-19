# Log-likelihood function for EPS-full
# No interactions, no confounding
# Hardy-Weinberg eq assumed

epsfullloglikhwe_maf_xe = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,n_cc,n_ic,
                                   genotypes,genoprobs,ng){

    len = 1 + dim(x_cc)[2] + dim(g_cc)[2]
    param = parameters

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp + exp(-(0.5*(y_ic - alpha - x_ic%*%betaE - c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs

    ind = rep(1,n_cc)
    for(j in 1:ng){
        ind = ind + g_cc[,j]*(3)^(j-1)
    }
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) - sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)

}

epsfullloglikhwe_maf = function(parameters,y_cc,y_ic,g_cc,n_cc,n_ic,
                                   genotypes,genoprobs,ng){

    len = 1 + dim(g_cc)[2]
    param = parameters

    lenp = length(param)
    alpha = param[1]
    betaG = param[2:len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp + exp(-(0.5*(y_ic - alpha -c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs
    ind = rep(1,n_cc)
    for(j in 1:ng){
        ind = ind + g_cc[,j]*(3)^(j-1)
    }
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) - sum(0.5*(y_cc - alpha - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)

}

epsfullloglikhwe_xe = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,n_cc,n_ic,ng){

    len = 1 + dim(x_cc)[2] + dim(g_cc)[2]

    m = parameters[(length(parameters)-ng + 1):length(parameters)]
    maf = exp(m)/(1+exp(m))
    param = parameters[1:(length(parameters)-ng)]

    getgeno = genprobhwe(ng,maf,geneffect = "additive")
    genotypes = as.matrix(getgeno[[1]])
    genoprobs = getgeno[[2]]

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp + exp(-(0.5*(y_ic - alpha - x_ic%*%betaE - c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs

    ind = rep(1,n_cc)
    for(j in 1:ng){
        ind = ind + g_cc[,j]*(3)^(j-1)
    }
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) - sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)

}

epsfullloglikhwe = function(parameters,y_cc,y_ic,g_cc,n_cc,n_ic,ng){

    len = 1 + dim(g_cc)[2]

    m = parameters[(length(parameters)-ng + 1):length(parameters)]
    maf = exp(m)/(1+exp(m))
    param = parameters[1:(length(parameters)-ng)]

    getgeno = genprobhwe(ng,maf,geneffect = "additive")
    genotypes = as.matrix(getgeno[[1]])
    genoprobs = getgeno[[2]]

    lenp = length(param)
    alpha = param[1]
    betaG = param[2:len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (i in 1:dim(genotypes)[1]){
        temp = temp + exp(-(0.5*(y_ic - alpha -c(t(genotypes[i,])%*%betaG))^2/sigma2))*genoprobs[i]
    }
    temp = sum(log(temp*(1/(sqrt(2*pi)*sigma))))

    # Add up contribution of completely observed individuals
    # with respect to the log of probability of their individual genotypes.
    # Probabilities extracted from genoprobs
    ind = rep(1,n_cc)
    for(j in 1:ng){
        ind = ind + g_cc[,j]*(3)^(j-1)
    }
    temp2 = sum(log(genoprobs[ind]))

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) - sum(0.5*(y_cc - alpha - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)

}
