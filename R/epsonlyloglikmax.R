# Maximize log-likelihood function for EPS-only data

epsonlyloglikmax = function(data,cutoffs,ll= FALSE, hessian = FALSE ){
    if(dim(data)[2]>1){
        y = data[,1]
        len = dim(data)[2]
        x = as.matrix(data[,2:len])
        init = lm(y ~ x)
        a.init = c(init$coef,log(summary(init)$sigma^2))
        result = optim(a.init,
                       epsonlyloglik,
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
                       epsonlyloglik_y,
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

}
