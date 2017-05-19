# Maximimize log-likelihood function for EPS-full
# No interactions, no confounding

epsfullloglikmax = function(data,ng,hwe = FALSE,maf=NA,
                            ll= FALSE, hessian = FALSE){
    if(hwe){
        ###############################################################
        # Hardy-Weinberg assumed, independent SNPs
        ###############################################################
        if(!is.na(as.matrix(maf)[1,1])){
            ###############################################################
            # MAF given by user
            ###############################################################
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

                getgeno = genprobhwe(ng,maf,geneffect = "additive")
                genotypes = as.matrix(getgeno[[1]])
                genoprobs = getgeno[[2]]

                init = lm(y_cc ~ x_cc + g_cc)

                a.init = c(init$coef,log(summary(init)$sigma^2))

                result = optim(a.init,
                               epsfullloglikhwe_maf_xe,
                               y_cc = y_cc,
                               y_ic = y_ic,
                               x_cc = x_cc,
                               x_ic = x_ic,
                               g_cc = g_cc,
                               n_cc = n_cc,
                               n_ic = n_ic,
                               genotypes = genotypes,
                               genoprobs = genoprobs,
                               ng = ng,
                               method = "BFGS",
                               control = list(fnscale = -1,
                                              trace = 1),
                               hessian = TRUE)
                mles = result$par[1:(length(result$par) - 1)]
                nparam = length(mles)
                sigma = sqrt(exp(result$par[(nparam+1)]))
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

                getgeno = genprobhwe(ng,maf,geneffect = "additive")
                genotypes = as.matrix(getgeno[[1]])
                genoprobs = getgeno[[2]]

                init = lm(y_cc ~ g_cc)

                a.init = c(init$coef,log(summary(init)$sigma^2))

                result = optim(a.init,
                               epsfullloglikhwe_maf,
                               y_cc = y_cc,
                               y_ic = y_ic,
                               g_cc = g_cc,
                               n_cc = n_cc,
                               n_ic = n_ic,
                               genotypes = genotypes,
                               genoprobs = genoprobs,
                               ng = ng,
                               method = "BFGS",
                               control = list(fnscale = -1,
                                              trace = 1),
                               hessian = TRUE)
                mles = result$par[1:(length(result$par) - 1)]
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
        }else{
            ###############################################################
            # MAF unknown
            ###############################################################
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

                init = lm(y_cc ~ x_cc + g_cc)

                mafinit = rep(0.2,ng)
                a.init = c(init$coef,log(summary(init)$sigma^2),log(mafinit/(1-mafinit)))

                result = optim(a.init,
                               epsfullloglikhwe_xe,
                               y_cc = y_cc,
                               y_ic = y_ic,
                               x_cc = x_cc,
                               x_ic = x_ic,
                               g_cc = g_cc,
                               n_cc = n_cc,
                               n_ic = n_ic,
                               ng = ng,
                               method = "BFGS",
                               control = list(fnscale = -1,
                                              trace = 1),
                               hessian = TRUE)
                mles = result$par[1:(length(result$par) -(2*ng) - 1)]
                nparam = length(mles)
                sigma = sqrt(exp(result$par[(nparam+1)]))
                m = result$par[(nparam+2):(length(result$par))]
                maf = exp(m)/(1+exp(m))
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

                init = lm(y_cc ~ g_cc)

                mafinit = rep(0.2,ng)
                a.init = c(init$coef,log(summary(init)$sigma^2),log(mafinit/(1-mafinit)))

                result = optim(a.init,
                               epsfullloglikhwe,
                               y_cc = y_cc,
                               y_ic = y_ic,
                               g_cc = g_cc,
                               n_cc = n_cc,
                               n_ic = n_ic,
                               ng = ng,
                               method = "BFGS",
                               control = list(fnscale = -1,
                                              trace = 1),
                               hessian = TRUE)
                mles = result$par[1:(length(result$par) -(2*ng) - 1)]
                nparam = length(mles)
                sigma = sqrt(exp(result$par[(nparam+1)]))
                m = result$par[(nparam+2):(length(result$par))]
                maf = exp(m)/(1+exp(m))
            }
            if(ll){
                return(list(result$value, c(mles,sigma,maf)))
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
                           epsfullloglik_xe,
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
                           epsfullloglik,
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
}
