# Log-likelihood function for EPS-full
# No interactions, no confounding
# Genotype frequencies known

epsfullloglikgfreq = function(parameters,data,ng,gfreq,geneffect = "additive"){
    len = dim(data)[2]
    param = parameters[1:(len+1)]

    # Extract all combinations of genotypes, and the
    # corresponding probability for each combination.
    # It is assumed that all genetic variables are statistically
    # independent. Use function genprob().

    temp1 = list()
    temp2 = list()

    index = 1
    for(i in 1:ng){
        temp1[[i]] = c(0,1,2)
        temp2[[i]] = c(gfreq[index],gfreq[index+1],gfreq[index+2])
        index = index + 3
    }
    geno = expand.grid(temp1)
    probs = expand.grid(temp2)
    probs2 = apply(probs,1,prod)

    genotypes = as.matrix(geno)
    genoprobs = probs2

    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]

    if(len > (ng + 1)){
        ###############################################################
        # Environmental covariates (xe) present
        ###############################################################
        y_cc = data_cc[,1]
        y_ic = data_ic[,1]
        g_cc = as.matrix(data_cc[,(len-ng+1):len])
        x_cc = as.matrix(data_cc[,2:(len-ng)])
        x_ic = as.matrix(data_ic[,2:(len-ng)])
        n_cc = length(y_cc)
        n_ic = length(y_ic)

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
        if(geneffect=="additive"){rowind = 3}else{rowind = 2}
        ind = rep(1,n_cc)
        for(j in 1:ng){
            ind = ind + g_cc[,j]*(rowind)^(j-1)
        }
        temp2 = sum(log(genoprobs[ind]))

        # Add up contribution of completely observed individuals
        # with respect to the conditional distribution of y
        temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
            sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

        # The log-likelihood
        ll = temp2 + temp3 + temp
        return(ll)
    }else{
        ###############################################################
        # Environmental covariates (xe) not present
        ###############################################################
        y_cc = data_cc[,1]
        y_ic = data_ic
        g_cc = as.matrix(data_cc[,(len-ng+1):len])
        n_cc = length(y_cc)
        n_ic = length(y_ic)

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
            sum(0.5*(y_cc - alpha - g_cc%*%betaG)^2/sigma2)

        # The log-likelihood
        ll = temp2 + temp3 + temp
        return(ll)
    }
}
