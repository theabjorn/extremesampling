# Log-likelihood function for EPS-only

epsonlyloglik = function(parameters,data,cutoffs){
    param = parameters
    l = min(cutoffs)
    u = max(cutoffs)
    len = dim(data)[2]

    y = data[,1]
    x = as.matrix(data[,2:len])
    n = length(y)

    lenp = length(param)
    alpha = param[1]
    beta = param[2:(lenp-1)]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    z = (y-alpha - x%*%beta)/sigma
    zl = (l-alpha - x%*%beta)/sigma
    zu = (u-alpha - x%*%beta)/sigma
    ll = sum(log(dnorm(z)/sigma)-log(1-pnorm(zu)+pnorm(zl)))
    return(ll)
}

epsonlyloglik_y = function(parameters,data,cutoffs){
    param = parameters
    l = min(cutoffs)
    u = max(cutoffs)
    len = dim(data)[2]

    y = data[,1]
    n = length(y)

    lenp = length(param)
    alpha = param[1]
    tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)

    z = (y-alpha )/sigma
    zl = (l-alpha )/sigma
    zu = (u-alpha )/sigma
    ll = sum(log(dnorm(z)/sigma)-log(1-pnorm(zu)+pnorm(zl)))
    return(ll)
}
