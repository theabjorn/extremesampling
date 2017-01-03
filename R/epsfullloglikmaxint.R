# Maximimize log-likelihood function for EPS-full
# Interactions, no confounding

epsfullloglikmaxint = function(data,ng,interactind = list(NA),hwe = FALSE,
                               maf=NA,geneffect="additive",
                               ll= FALSE, hessian = FALSE){
    y = data[,1]
    xe = as.matrix(data[,2:(dim(data)[2]-ng)])
    xg = as.matrix(data[,(dim(data)[2]-ng+1):(dim(data)[2])])
    if(is.na(interactind[[1]][1])){
        interactind = list()
        for(i in 1:dim(xg)[2]){
            for(j in 1:dim(xe)[2]){
                interactind[[(length(interactind)+1)]] = c(i,j)
            }
        }
    }
    xeg = matrix(NA,ncol = length(interactind),nrow = dim(xe)[1])

    for(i in 1:length(interactind)){
        xeg[,i] = xe[,interactind[[i]][2]]*xg[,interactind[[i]][1]]
    }
    alldata = cbind(data,xeg)

    if(hwe){
        ###############################################################
        # Hardy-Weinberg assumed
        ###############################################################
        if(!is.na(as.matrix(maf)[1,1])){
            ###############################################################
            # MAF given by user
            ###############################################################
            len = dim(alldata)[2]
            data_cc = alldata[!is.na(alldata[,len]),]
            y = data_cc[,1]
            x = as.matrix(data_cc[,2:len])
            init = lm(y ~ x)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsfullloglikinthwe,
                           data = data,
                           ng = ng,
                           interactind = interactind,
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
            len = dim(alldata)[2]
            data_cc = alldata[!is.na(alldata[,len]),]
            y = data_cc[,1]
            x = as.matrix(data_cc[,2:len])
            init = lm(y ~ x)
            mafinit = rep(0.2,ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),
                       log(mafinit/(1-mafinit)))
            result = optim(a.init,
                           epsfullloglikinthwe,
                           data = data,
                           ng = ng,
                           interactind = interactind,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)

            mles = result$par[1:(length(result$par)-ng - 1)]
            sigma = sqrt(exp(result$par[(length(result$par)-ng)]))
            ind = length(result$par)-ng
            mafs = exp(result$par[(ind + 1):(length(result$par))])/
                (1+exp(result$par[(ind + 1):(length(result$par))]))
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
        len = dim(alldata)[2]
        data_cc = alldata[!is.na(alldata[,len]),]
        y = data_cc[,1]
        x = as.matrix(data_cc[,2:len])
        init = lm(y ~ x)
        if(geneffect == "additive"){
#             count = 1
#             probinit = c()
#             for(c in 1:ng){
#                 probinit[count] = log(0.2/(1-0.2))
#                 probinit[count+1] = log(0.2/(1-probinit[count]-0.2))
#             }
            probinit = rep(0.2,2*ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),
                       log(probinit/(1-probinit)))
            result = optim(a.init,
                           epsfullloglikint,
                           data = data,
                           ng = ng,
                           interactind = interactind,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)-2*ng - 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))
            probabilities = c()
            count = nparam + 2
            count2 = 1
            for(r in 1:(ng)){
                probabilities[count2] = exp(result$par[count])/
                    (1+exp(result$par[count]))
                probabilities[count2+1] = (1-probabilities[count2])*
                    exp(result$par[count+1])/(1+exp(result$par[count+1]))
                probabilities[count2+2] = 1 - probabilities[count2] -
                    probabilities[count2+1]
                count2 = count2 + 3
                count = count + 2
            }
        }else{
            probinit = rep(0.2,ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),
                       log(probinit/(1-probinit)))
            result = optim(a.init,
                           epsfullloglikint,
                           data = data,
                           ng = ng,
                           interactind = interactind,
                           geneffect = geneffect,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)-ng - 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))
            probabilities = c()
            count = nparam + 2
            count2 = 1
            for(r in 1:(ng)){
                probabilities[count2] = exp(result$par[count])/
                    (1+exp(result$par[count]))
                probabilities[count2+1] = 1 - exp(result$par[count])/
                    (1+exp(result$par[count]))
                count2 = count2 + 2
                count = count + 1
            }
        }
        if(ll){
            return(list(result$value, c(mles,sigma)))
        }else if(hessian){
            return(list(result$hessian,c(mles,sigma,probabilities)))
        }else{
            return(c(mles,sigma,probabilities))
        }
    }
}
