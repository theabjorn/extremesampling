# Maximimize log-likelihood function for EPS-complete
# No interactions, no confounding

epscomploglikmax = function(data,ng,hwe = FALSE,maf=NA,gfreq=NA,
                            geneffect="additive",
                            ll= FALSE, hessian = FALSE){
    if(hwe){
        ###############################################################
        # Hardy-Weinberg assumed
        ###############################################################
        if(!is.na(as.matrix(maf)[1,1])){
            ###############################################################
            # MAF given by user
            ###############################################################
            len = dim(data)[2]
            data_cc = data[!is.na(data[,len]),]
            y = data_cc[,1]
            x = as.matrix(data_cc[,2:len])
            init = lm(y ~ x)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epscomploglikhwe,
                           data = data,
                           ng = ng,
                           maf = maf,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)-1)]
            sigma = sqrt(exp(result$par[(length(result$par))]))
            if(ll){
                return(list(result$value, c(mles,sigma)))
            }else if(hessian){
                return(list(result$hessian,c(mles,sigma)))
            }else{
                return(c(mles,sigma))
            }
        }else{
            ###############################################################
            # MAF unknown
            ###############################################################
            len = dim(data)[2]
            data_cc = data[!is.na(data[,len]),]
            y = data_cc[,1]
            x = as.matrix(data_cc[,2:len])
            init = lm(y ~ x)
            mafinit = rep(0.2,ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),
                       log(mafinit/(1-mafinit)))
            result = optim(a.init,
                           epscomploglikhwe,
                           data = data,
                           ng = ng,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)

            mles = result$par[1:(length(result$par)-ng - 1)]
            sigma = sqrt(exp(result$par[(length(result$par)-ng)]))
            mafs = exp(result$par[(length(result$par)-ng + 1):
                                      (length(result$par))])/
                (1+exp(result$par[(length(result$par)-ng + 1):
                                      (length(result$par))]))
            if(ll){
                return(list(result$value, c(mles,sigma,mafs)))
            }else if(hessian){
                return(list(result$hessian,c(mles,sigma,mafs)))
            }else{
                return(c(mles,sigma,mafs))
            }
        }
    }else{
        ###############################################################
        # Hardy-Weinberg not assumed
        ###############################################################
        len = dim(data)[2]
        data_cc = data[!is.na(data[,len]),]
        y = data_cc[,1]
        x = as.matrix(data_cc[,2:len])
        init = lm(y ~ x)

        if(is.na(gfreq[1])){
            probinit = rep(0.2,2*ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),
                       log(probinit/(1-probinit)))
            result = optim(a.init,
                           epscomploglik,
                           data = data,
                           ng = ng,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)-2*ng - 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))
        }else{
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epscomploglikgfreq,
                           data = data,
                           ng = ng,
                           gfreq = gfreq,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)- 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))
        }

        if(ll){
            return(list(result$value, c(mles,sigma)))
        }else if(hessian){
            return(list(result$hessian,c(mles,sigma)))
        }else{
            return(c(mles,sigma))
        }
    }
}
