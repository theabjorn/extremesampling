# Log-likelihood function for EPS-CC

epsCC.loglik_ex = function(parameters,data,l,u){
    param = parameters
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

epsCC.loglik_e = function(parameters,data,l,u){
    param = parameters
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

# epsCC.loglik_rand = function(parameters,data,l,u,randomindex){
#     param = parameters
#     len = dim(data)[2]
#
#     y_e = data[randomindex==0,1]
#     x_e = as.matrix(data[randomindex==0,2:len])
#     ne = length(y_e)
#     y_r = data[,1][randomindex==1]
#     x_r = as.matrix(data[randomindex==1,2:len])
#     nr = length(y_r)
#
#     lenp = length(param)
#     alpha = param[1]
#     beta = param[2:(lenp-1)]
#     tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)
#
#     z = (y_e-alpha - x_e%*%beta)/sigma
#     zl = (l-alpha - x_e%*%beta)/sigma
#     zu = (u-alpha - x_e%*%beta)/sigma
#
#     z_r = (y_r-alpha - x_r%*%beta)/sigma
#
#     ll = sum(log(dnorm(z_r)/sigma)) + sum(log(dnorm(z)/sigma)-log(1-pnorm(zu)+pnorm(zl)))
#     return(ll)
# }


#
# epsCC.loglik_z_rand = function(parameters,data,cutoffs,randomindex){
#     param = parameters
#     l = min(cutoffs)
#     u = max(cutoffs)
#
#     y_e = data[,1][randomindex==0]
#     ne = length(y_e)
#     y_r = data[,1][randomindex==1]
#     nr = length(y_r)
#
#     lenp = length(param)
#     alpha = param[1]
#     tau = param[lenp]; sigma2 = exp(tau); sigma = sqrt(sigma2)
#
#     z = (y_e-alpha)/sigma
#     zl = (l-alpha)/sigma
#     zu = (u-alpha)/sigma
#
#     z_r = (y_r-alpha)/sigma
#
#     ll = sum(log(dnorm(z_r)/sigma)) + sum(log(dnorm(z)/sigma)-log(1-pnorm(zu)+pnorm(zl)))
#     return(ll)
# }
#
#
# #######################################################
# # EPS-CC loglik for secondary phenotype W
# #######################################################
#
# epsCC.loglik.W = function(parameters,data){
#     param = parameters
#     len = dim(data)[2]
#
#     w = data[,1]
#     x = as.matrix(data[,2:(len-1)])
#     gamma = data[,len]
#     n = length(w)
#
#     lenp = length(param)
#     alpha = param[1]
#     beta = param[2:(lenp-2)]
#     tau = param[(lenp)-1]; sigma2 = exp(tau); sigma = sqrt(sigma2)
#     tmp = param[(lenp)]; rho = exp(tmp)/(1+exp(tmp))
#
#     z = (w-alpha - x%*%beta - sigma*rho*gamma)
#
#     ll = -(n/2)*log(sigma2) - (n/2)*log(1-rho^2) - 0.5*(1/sigma2)*(1/(1-rho^2))*sum(z^2)
#     return(ll)
# }
