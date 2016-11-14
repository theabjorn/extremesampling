# Maximimize log-likelihood function for EPS-complete
# Interactions, confounding


epscomploglikmaxcondint = function(data, ng, cind,interactind,
                                ll = FALSE, hessian = FALSE){

    y = data[,1]
    xe = as.matrix(data[,2:(dim(data)[2]-ng)])
    xg = as.matrix(data[,(dim(data)[2]-ng+1):(dim(data)[2])])
    xeg = matrix(NA,ncol = length(interactind),nrow = dim(xe)[1])

    for(i in 1:length(interactind)){
        xeg[,i] = xe[,interactind[[i]][2]]*xg[,interactind[[i]][1]]
    }
    alldata = cbind(data,xeg)

    len = dim(alldata)[2]
    data_cc = alldata[!is.na(alldata[,len]),]
    y = data_cc[,1]
    x = as.matrix(data_cc[,2:len])
    init = lm(y ~ x)

    xe = as.matrix(data_cc[,2:(len-ng)])
    xg = as.matrix(data_cc[,(len-ng+1):len])
    xu = as.matrix(unique(xe[,cind]))
    nu = dim(xu)[1]
    # For each unique xe, two probabilities per genetic variant
    nprob = 2*nu*ng
    probinit = rep(0.2,nprob)
    a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))
    result = optim(a.init,
                   epscomploglikintcond,
                   data = data,
                   ng = ng,
                   interactind = interactind,
                   cind = cind,
                   method = "BFGS",
                   control = list(fnscale = -1,
                                  trace = 1),
                   hessian = TRUE)
    ne = dim(xe)[2]
    nparam = 1+ne+ng+length(interactind)
    mles = result$par[1:nparam]
    sigma = sqrt(exp(result$par[(nparam+1)]))
    probs = c()
    count = nparam + 2
    count2 = 1
    for(r in 1:(nu*ng)){
      probs[count2] = exp(result$par[count])/(1+exp(result$par[count]))
      probs[count2+1] = (1-probs[count2])*exp(result$par[count+1])/
          (1+exp(result$par[count+1]))
      probs[count2+2] = 1 - probs[count2] - probs[count2+1]
      count2 = count2 + 3
      count = count + 2
    }
    if(ll){
      return(list(result$value, c(mles,sigma)))
    }else if(hessian){
      return(list(result$hessian,c(mles,sigma)))
    }else{
      return(c(mles,sigma,probs))
    }
}
