# Log-likelihood function for EPS-full
# Interactions, no confounding
# Hardy-Weinberg eq assumed

epsfullloglikinthwe = function(parameters,data,ng,interactind,maf = NA,
                               geneffect = "additive"){
    if(!is.na(as.matrix(maf)[1,1])){
        # MAF supplied by the user
        param = parameters
    }else{
        # MAF unknown (part of parameter vector)
        m = parameters[(length(parameters)-ng + 1):length(parameters)]
        maf = exp(m)/(1+exp(m))
        param = parameters[1:(length(parameters)-ng)]
    }
    # Extract all combinations of genotypes, and the
    # corresponding probability for each combination.
    # It is assumed that all genetic variables are statistically
    # independent. Use function genprobhwe().
    getgeno = genprobhwe(ng,maf,geneffect)
    genotypes = as.matrix(getgeno[[1]])
    genoprobs = getgeno[[2]]

    len = dim(data)[2]
    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]

    y_cc = data_cc[,1]
    y_ic = data_ic[,1]
    g_cc = as.matrix(data_cc[,(len-ng+1):len])
    x_cc = as.matrix(data_cc[,2:(len-ng)])
    x_ic = as.matrix(data_ic[,2:(len-ng)])
    n_cc = length(y_cc)
    n_ic = length(y_ic)

    # Interactions for complete case individuals
    neg = length(interactind)
    xg_cc = matrix(NA,ncol = length(interactind),nrow = n_cc)
    for(i in 1:neg){
        xg_cc[,i] = x_cc[,interactind[[i]][2]]*g_cc[,interactind[[i]][1]]
    }

    # Interactions (all possible) for missing-genotype individuals
    # Each list element represents
    intgenotypes = list()
    for(g in 1:dim(genotypes)[1]){
        intgenotypes[[g]] = matrix(NA,ncol = length(interactind),
                                        nrow = n_ic)
        for(j in 1:neg){
            intgenotypes[[g]][,j] =
                x_ic[,interactind[[j]][2]]*genotypes[g,interactind[[j]][1]]
        }
    }

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
                           intgenotypes[[g]]%*%betaEG)^2)/sigma2)*genoprobs[g]
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
