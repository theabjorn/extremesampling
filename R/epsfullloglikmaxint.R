# Maximimize log-likelihood function for EPS-full
# Interactions, no confounding

epsfullloglikmaxint = function(data,ng,interactind,hwe = FALSE,
                               maf = NA, ll = FALSE, hessian = FALSE){

    len = dim(data)[2]

    neg = length(interactind)

    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]

    y_cc = data_cc[,1]
    y_ic = data_ic[,1]
    g_cc = as.matrix(data_cc[,(len-ng+1):len])
    x_cc = as.matrix(data_cc[,2:(len-ng)])
    x_ic = as.matrix(data_ic[,2:(len-ng)])
    n_cc = length(y_cc)
    n_ic = length(y_ic)

    genotypes = unique(g_cc)

    xg_cc = matrix(NA,ncol = neg,nrow = n_cc)
    for(i in 1:neg){
        xg_cc[,i] = x_cc[,interactind[[i]][2]]*g_cc[,interactind[[i]][1]]
    }
    interactgenotypes = list()
    for (i in 1:dim(genotypes)[1]){
        interactgenotypes[[i]] = matrix(NA,ncol = length(interactind),nrow = n_ic)
        for(j in 1:neg){
            interactgenotypes[[i]][,j] = x_ic[,interactind[[j]][2]]*genotypes[i,interactind[[j]][1]]
        }
    }

    init = lm(y_cc ~ x_cc+g_cc+xg_cc)

    if(hwe){
        ###############################################################
        # Hardy-Weinberg assumed, independent SNPs
        ###############################################################
        if(!is.na(as.matrix(maf)[1,1])){
            ###############################################################
            # MAF given by user
            ###############################################################

            getgeno = genprobhwe(ng,maf,geneffect = "additive")
            genotypes = as.matrix(getgeno[[1]])
            genoprobs = getgeno[[2]]

            a.init = c(init$coef,log(summary(init)$sigma^2))

            result = optim(a.init,
                           epsfullloglikinthwe_maf,
                           y_cc = y_cc,
                           y_ic = y_ic,
                           x_cc = x_cc,
                           x_ic = x_ic,
                           g_cc = g_cc,
                           xg_cc = xg_cc,
                           n_cc = n_cc,
                           n_ic = n_ic,
                           interactind = interactind,
                           interactgenotypes = interactgenotypes,
                           maf = maf,
                           genotypes = genotypes,
                           genoprobs = genoprobs,
                           ng = ng,
                           neg = neg,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par) - 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))

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
            probinit = rep(0.2,ng)
            a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))

            result = optim(a.init,
                           epsfullloglikinthwe,
                           y_cc = y_cc,
                           y_ic = y_ic,
                           x_cc = x_cc,
                           x_ic = x_ic,
                           g_cc = g_cc,
                           xg_cc = xg_cc,
                           n_cc = n_cc,
                           n_ic = n_ic,
                           interactind = interactind,
                           interactgenotypes = interactgenotypes,
                           ng = ng,
                           neg = neg,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            mles = result$par[1:(length(result$par)-(2*ng) - 1)]
            nparam = length(mles)
            sigma = sqrt(exp(result$par[(nparam+1)]))
            m = result$par[(nparam+2):(length(result$par))]
            maf = exp(m)/(1+exp(m))

            if(ll){
                return(list(result$value, c(mles,sigma)))
            }else if(hessian){
                return(list(result$hessian,c(mles,sigma,maf)))
            }else{
                return(c(mles,sigma,maf))
            }
        }
    }else{
        ###############################################################
        # Hardy-Weinberg not assumed
        ###############################################################
        nug = dim(unique(g_cc))[1]

        probinit = rep(1/(nug+1),(nug-1))
        a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))

        result = optim(a.init,
                       epsfullloglikint,
                       y_cc = y_cc,
                       y_ic = y_ic,
                       x_cc = x_cc,
                       x_ic = x_ic,
                       g_cc = g_cc,
                       xg_cc = xg_cc,
                       n_cc = n_cc,
                       n_ic = n_ic,
                       genotypes = genotypes,
                       interactind = interactind,
                       interactgenotypes = interactgenotypes,
                       ng = ng,
                       neg = neg,
                       method = "BFGS",
                       control = list(fnscale = -1,
                                      trace = 1),
                       hessian = TRUE)
        mles = result$par[1:(length(result$par)-(nug-1) - 1)]
        nparam = length(mles)
        sigma = sqrt(exp(result$par[(nparam+1)]))
        m = result$par[(nparam+2):(length(result$par))]
        probs = c()
        probs[1] = exp(m[1])/(1+exp(m[1]))
        for(k in 2:(nug-1)){
            probs[k] = (1-sum(probs))*exp(m[k])/(1+exp(m[k]))
        }
        gfreqs = c(probs,1-sum(probs))

        if(ll){
            return(list(result$value, c(mles,sigma),gfreqs))
        }else if(hessian){
            return(list(result$hessian,c(mles,sigma),gfreqs,genotypes))
        }else{
            return(c(mles,sigma,gfreqs))
        }
    }
}
