# Log-likelihood function for EPS-complete
# Confounding assumed, no interactions

epscomploglikcond = function(parameters,data,ng,cind,geneffect = "additive"){
    len = dim(data)[2]
    m = parameters[(len+2):length(parameters)]
    param = parameters[1:(len+1)]

    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]
    y_cc = data_cc[,1]
    y_ic = data_ic[,1]
    g_cc = as.matrix(data_cc[,(len-ng+1):len])
    x_cc = as.matrix(data_cc[,2:(len-ng)])
    x_ic = as.matrix(data_ic[,2:(len-ng)])
    n_cc = length(y_cc)
    n_ic = length(y_ic)

    xu = as.matrix(unique(x_cc[,cind]))
    nu = dim(xu)[1]

    allgenotypes = list()
    allgenoprobs = list()

    # For each unique level of xe, each genetic covariate has a
    # probability distribution, so that each unique combination of
    # genotypes has a probability. These are stored in allgenotypes and
    # allgenoprobs

    if(geneffect == "additive"){
        m2 = m
        for(i in 1:nu){
            # For each level of xe, get probabilities
            # for each combination of genotypes
            probs = c()
            count = 1
            for(t in 1:ng){
                probs[count] = exp(m2[count])/(1+exp(m2[count]))
                probs[count + 1] = (1-probs[count])*
                    exp(m2[count+1])/(1+exp(m2[count+1]))
                count = count + 2
            }
            getgeno = genprob(ng,probs)
            allgenotypes[[i]] = as.matrix(getgeno[[1]])
            allgenoprobs[[i]] = getgeno[[2]]
            m2 = m2[count:length(m2)]
        }
    }else{
        m2 = m
        for(i in 1:nu){
            probs = exp(m2[1:ng])/(1+exp(m2[1:ng]))
            getgeno = genprob(ng,probs)
            allgenotypes[[i]] = as.matrix(getgeno[[1]])
            allgenoprobs[[i]] = getgeno[[2]]
            m2 = m2[(ng+1):length(m2)]
        }
    }

    #print(allgenoprobs)

    lenp = length(param)
    alpha = param[1]
    betaE = param[2:(len-ng)]
    betaG = param[(len-ng+1):len]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    temp = 0
    temp2 = 0
    if(geneffect=="additive"){rowind = 3}else{rowind = 2}

    for(k in 1:nu){
        # Add up contribution to log-likelihood for all individuals
        # with missing genotypes, summing over all combinations of genotypes
        # Blocks of individuals with the same level of xe are
        # considered simultaneously
        idrows = colSums((data.frame(t(x_ic[,cind])) == t(xu[k,]))) # OBS x_ic
        id = idrows > 0 # index for individuals in block k
        genotypes = allgenotypes[[k]]
        genoprobs = allgenoprobs[[k]]
        tempk = 0
        for (j in 1:dim(genotypes)[1]){
            tempk = tempk +
                exp(-(0.5*(y_ic[id] - alpha -
                               as.matrix(x_ic[id,])%*%betaE -
                               c(t(genotypes[j,])%*%betaG))^2/sigma2)
                    )*genoprobs[j]
        }
        tempk = sum(log(tempk*(1/(sqrt(2*pi)*sigma))))
        temp = temp + tempk

        # Add up contribution of completely observed individuals
        # with respect to the log of probability of their individual genotypes.
        # Blocks of individuals with the same level of xe are
        # considered simultaneously
        # Probabilities extracted from genoprobs
        idrows = colSums((data.frame(t(x_cc[,cind])) == t(xu[k,]))) # OBS x_cc
        id = idrows > 0 # index for individuals in block k
        genotypes = allgenotypes[[k]]
        genoprobs = allgenoprobs[[k]]
        n_cck = sum(id)
        ind = rep(1,n_cck)
        for(j in 1:ng){
            ind = ind + g_cc[id,j]*(rowind)^(j-1)
        }
        temp2 = temp2 + sum(log(genoprobs[ind]))

    }

    # Add up contribution of completely observed individuals
    # with respect to the conditional distribution of y
    temp3 = - 0.5*n_cc*log(2*pi) - n_cc*log(sigma) -
        sum(0.5*(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2/sigma2)

    # The log-likelihood
    ll = temp2 + temp3 + temp
    return(ll)
}
