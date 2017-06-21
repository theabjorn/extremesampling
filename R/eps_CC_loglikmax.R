# Maximize log-likelihood function for EPS-only data

epsCC.loglikmax = function(data,cutoffs,randomindex,ll=FALSE, hessian = FALSE){

    rsample = TRUE
    if(missing(randomindex)){
        rsample = FALSE
    }else if(sum(randomindex)==0){
        rsample = FALSE
    }

    if(!rsample){
        if(dim(data)[2]>1){
            y = data[,1]
            len = dim(data)[2]
            x = as.matrix(data[,2:len])
            init = lm(y ~ x)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik,
                           data = data,
                           cutoffs = cutoffs,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }else{
            y = data[,1]
            init = lm(y ~ 1)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik_z,
                           data = data,
                           cutoffs = cutoffs,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }

    }else{
        if(dim(data)[2]>1){
            y = data[,1]
            len = dim(data)[2]
            x = as.matrix(data[,2:len])
            init = lm(y ~ x)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik_rand,
                           data = data,
                           cutoffs = cutoffs,
                           randomindex = randomindex,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }else{
            y = data[,1]
            init = lm(y ~ 1)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik_z_rand,
                           data = data,
                           cutoffs = cutoffs,
                           randomindex=randomindex,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }
    }



}
