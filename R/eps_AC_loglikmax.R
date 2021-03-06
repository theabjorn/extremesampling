# Maximimize log-likelihood function for EPS-full
# No confounding

epsAC.loglikmax = function(data,ng,ll= FALSE, hessian = FALSE){

        if(dim(data)[2] > (1+ng)){
            ###############################################################
            # Environmental covariates (xe) present
            ###############################################################
            len = dim(data)[2]

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
            nug = dim(genotypes)[1]

            init = lm(y_cc ~ x_cc + g_cc)

            probinit = rep(1/(2*nug),(nug-1))
            a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))

            result = optim(a.init,
                           epsAC.loglik_xe,
                           y_cc = y_cc,
                           y_ic = y_ic,
                           x_cc = x_cc,
                           x_ic = x_ic,
                           g_cc = g_cc,
                           n_cc = n_cc,
                           n_ic = n_ic,
                           genotypes = genotypes,
                           ng = ng,
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

        }else{
            ###############################################################
            # Environmental covariates not present
            ###############################################################
            len = dim(data)[2]

            data_cc = data[!is.na(data[,len]),]
            data_ic = data[is.na(data[,len]),1:(len-ng)]

            y_cc = data_cc[,1]
            y_ic = data_ic
            g_cc = as.matrix(data_cc[,(len-ng+1):len])
            n_cc = length(y_cc)
            n_ic = length(y_ic)

            genotypes = unique(g_cc)
            nug = dim(genotypes)[1]

            init = lm(y_cc ~ g_cc)

            probinit = rep(1/(2*nug),(nug-1))
            a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))

            result = optim(a.init,
                           epsAC.loglik,
                           y_cc = y_cc,
                           y_ic = y_ic,
                           g_cc = g_cc,
                           n_cc = n_cc,
                           n_ic = n_ic,
                           genotypes = genotypes,
                           ng = ng,
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
        }

        if(ll){
            return(list(result$value, c(mles,sigma),gfreqs))
        }else if(hessian){
            return(list(result$hessian,c(mles,sigma),gfreqs,genotypes))
        }else{
            return(c(mles,sigma,gfreqs))
        }
}

# Maximimize log-likelihood function for EPS-AC
# Confounding

epsAC.loglikmaxcond = function(data, confounder, ng, ll = FALSE, hessian = FALSE,snpnames){

    if(missing(snpnames)){
        gs = 1:ng
        snpnames = paste0("snp",gs)
    }

    len = dim(data)[2]

    data_cc = data[!is.na(data[,len]),]
    data_ic = data[is.na(data[,len]),1:(len-ng)]

    y_cc = data_cc[,1]
    y_ic = data_ic[,1]
    g_cc = as.matrix(data_cc[,(len-ng+1):len])
    x_cc = as.matrix(data_cc[,2:(len-ng)])
    x_ic = as.matrix(data_ic[,2:(len-ng)])
    n_cc = length(y_cc)
    n_ic = length(y_ic)
    confounder_ic = as.matrix(confounder)[is.na(data[,len]),]
    confounder_cc = as.matrix(confounder)[!is.na(data[,len]),]

    genotypes = unique(g_cc)
    nug = dim(genotypes)[1]

    init = lm(y_cc ~ x_cc + g_cc)

    xu = as.matrix(unique(confounder))
    nu = dim(xu)[1]

    # For each unique xe, two probabilities per genetic variant
    nprob = nu*(nug-1)

    probinit = rep(1/(2*nprob),nprob)
    a.init = c(init$coef,log(summary(init)$sigma^2),log(probinit/(1-probinit)))
    result = optim(a.init,
                   epsAC.loglikcond,
                   y_cc = y_cc,
                   y_ic = y_ic,
                   x_cc = x_cc,
                   x_ic = x_ic,
                   g_cc = g_cc,
                   n_cc = n_cc,
                   n_ic = n_ic,
                   genotypes = genotypes,
                   ng = ng,
                   xu = xu,
                   nu = nu,
                   confounder_cc = confounder_cc,
                   confounder_ic = confounder_ic,
                   method = "BFGS",
                   control = list(fnscale = -1,
                                  trace = 1),
                   hessian = TRUE)

    mles = result$par[1:(length(result$par)-(nprob) - 1)]
    nparam = length(mles)
    sigma = sqrt(exp(result$par[(nparam+1)]))
    m = result$par[(nparam+2):(length(result$par))]

    allgenoprobs = list()
    for(i in 1:nu){
        mtemp = m[1:(dim(genotypes)[1]-1)]
        m = m[dim(genotypes)[1]:length(m)]
        nug = length(mtemp)+1
        probs = c()
        probs[1] = exp(mtemp[1])/(1+exp(mtemp[1]))
        for(k in 2:(nug-1)){
            probs[k] = (1-sum(probs))*exp(mtemp[k])/(1+exp(mtemp[k]))
        }
        allgenoprobs[[i]] = cbind(c(probs,1-sum(probs)),genotypes)
        colnames(allgenoprobs[[i]]) = c(paste0("P(Xg|Xe=",paste(xu[i,],collapse = "-"),")"),snpnames)
    }
    if(ll){
        return(list(result$value, c(mles,sigma),allgenoprobs))
    }else if(hessian){
        return(list(result$hessian,c(mles,sigma),allgenoprobs))
    }else{
        return(c(mles,sigma))
    }
}









