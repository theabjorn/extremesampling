# Log-likelihood function for EPS-full
# Interactions, no confounding
# Hardy-Weinberg eq assumed

epsfullloglikinthwe_maf = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,xg_cc,n_cc,n_ic,
                                   interactind,interactgenotypes,maf,genotypes,genoprobs,ng,neg){

    len = 1 + dim(x_cc)[2] + dim(g_cc)[2]
    param = parameters

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    betaEG = param[(len+1):(lenp-1)]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (g in 1:dim(genotypes)[1]){
        temp = temp +
            exp(-(0.5*(y_ic - alpha - x_ic%*%betaE -
                           c(t(genotypes[g,])%*%betaG) -
                           interactgenotypes[[g]]%*%betaEG)^2)/sigma2)*genoprobs[g]
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
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - x_cc%*%betaE -
                     g_cc%*%betaG - xg_cc%*%betaEG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}



epsfullloglikinthwe = function(parameters,y_cc,y_ic,x_cc,x_ic,g_cc,xg_cc,n_cc,n_ic,
                               interactind,interactgenotypes,ng,neg){

    len = 1 + dim(x_cc)[2] + dim(g_cc)[2]

    m = parameters[(length(parameters)-ng + 1):length(parameters)]
    maf = exp(m)/(1+exp(m))
    param = parameters[1:(length(parameters)-ng)]

    # Extract all combinations of genotypes, and the
    # corresponding probability for each combination.
    # It is assumed that all genetic variables are statistically
    # independent. Use function genprobhwe().
    getgeno = genprobhwe(ng,maf,geneffect="additive")
    genotypes = as.matrix(getgeno[[1]])
    genoprobs = getgeno[[2]]

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    betaEG = param[(len+1):(lenp-1)]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    # Add up contribution to log-likelihood for all individuals
    # with missing genotypes, summing over all combinations of genotypes
    temp = 0
    for (g in 1:dim(genotypes)[1]){
        temp = temp +
            exp(-(0.5*(y_ic - alpha - x_ic%*%betaE -
                           c(t(genotypes[g,])%*%betaG) -
                           interactgenotypes[[g]]%*%betaEG)^2)/sigma2)*genoprobs[g]
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
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - x_cc%*%betaE -
                     g_cc%*%betaG - xg_cc%*%betaEG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}
