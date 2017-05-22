hfun0 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs){
    tmp = rep(0,length(y))
    tmp2 = rep(0,length(y))
    for (i in 1:dim(genotypes)[1]){
        fi = y - alpha - xe%*%betaE - c(t(genotypes[i,])%*%betaG)
        tmp = tmp + (fi^a)*dnorm(fi/sigma)*genoprobs[i]*(genotypes[i,gint]^c)
        tmp2 = tmp2 + dnorm(fi/sigma)*genoprobs[i]
    }
    return(tmp/tmp2)
}

hfun0dash = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,k){
    fk = y - alpha - xe%*%betaE - c(t(genotypes[k,])%*%betaG)
    mk = length(genoprobs)
    fmk = y - alpha - xe%*%betaE - c(t(genotypes[mk,])%*%betaG)

    tmp = (fk^a)*dnorm(fk/sigma)*(genotypes[k,gint]^c) - (fmk^a)*dnorm(fmk/sigma)*(genotypes[mk,gint]^c)

    tmp2 = rep(0,length(y))
    for (j in 1:dim(genotypes)[1]){
        fi = y - alpha - xe%*%betaE - c(t(genotypes[j,])%*%betaG)
        tmp2 = tmp2 + dnorm(fi/sigma)*genoprobs[j]
    }
    return(tmp/tmp2)
}

hfun1dash = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,k){
    fk = y - alpha - xe%*%betaE - c(t(genotypes[k,])%*%betaG)
    mk = length(genoprobs)
    fmk = y - alpha - xe%*%betaE - c(t(genotypes[mk,])%*%betaG)

    tmp = (fk^a)*dnorm(fk/sigma)*(genotypes[k,gint]^c)*genotypes[k,] - (fmk^a)*dnorm(fmk/sigma)*(genotypes[mk,gint]^c)*genotypes[mk,]

    tmp2 = rep(0,length(y))
    for (i in 1:dim(genotypes)[1]){
        fi = y - alpha - xe%*%betaE - c(t(genotypes[i,])%*%betaG)
        tmp2 = tmp2 + dnorm(fi/sigma)*genoprobs[i]
    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}

hfun1 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs){
    tmp = matrix(0,ncol=dim(genotypes)[2],nrow = length(y))
    tmp2 = rep(0,length(y))
    for (j in 1:dim(genotypes)[1]){
        fi = y - alpha - xe%*%betaE - c(t(genotypes[j,])%*%betaG)
        #tmp = tmp + ((fi^a)*dnorm(fi/sigma)*genoprobs[j]*(genotypes[j,gint]^c))*c(genotypes[j,])
        for(k in 1:dim(genotypes)[2]){
            tmp[,k] = tmp[,k] + ((fi^a)*dnorm(fi/sigma)*genoprobs[j]*(genotypes[j,gint]^c))*c(genotypes[j,k])
        }

        tmp2 = tmp2 + dnorm(fi/sigma)*genoprobs[j]
    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}

hfun2 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs){
    tmp = matrix(0,ncol=dim(genotypes)[2]*dim(genotypes)[2],nrow = length(y))
    tmp2 = rep(0,length(y))
    for (j in 1:dim(genotypes)[1]){
        fi = y - alpha - xe%*%betaE - c(t(genotypes[j,])%*%betaG)
        for(k in 1:(dim(genotypes)[2]*dim(genotypes)[2])){
            tmp[,k] = tmp[,k] + ((fi^a)*dnorm(fi/sigma)*genoprobs[j]*(genotypes[j,gint]^c))*(as.vector(t(c(genotypes[j,]%*%t(genotypes[j,]))))[k])
        }
        tmp2 = tmp2 + dnorm(fi/sigma)*genoprobs[j]
    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}


chfun0 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,uindex){
    nu = length(uindex)
    tmp = rep(0,length(y))
    tmp2 = rep(0,length(y))
    for(u in 1:nu){
        cy = y[(uindex[[u]]==1)]
        cxe = as.matrix(xe[(uindex[[u]]==1),])
        cgenoprobs = genoprobs[[u]][,1]
        for (i in 1:dim(genotypes)[1]){
            fi = cy - alpha - cxe%*%betaE - c(t(genotypes[i,])%*%betaG)
            tmp[(uindex[[u]]==1)] = tmp[(uindex[[u]]==1)] + (fi^a)*dnorm(fi/sigma)*cgenoprobs[i]*(genotypes[i,gint]^c)
            tmp2[(uindex[[u]]==1)] = tmp2[(uindex[[u]]==1)] + dnorm(fi/sigma)*cgenoprobs[i]
        }
    }
    return(tmp/tmp2)
}

chfun0dash = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,k,uindex,u1){
    nu = length(uindex)

    fk = y - alpha - xe%*%betaE - c(t(genotypes[k,])%*%betaG)
    mk = dim(genotypes)[1]
    fmk = y - alpha - xe%*%betaE - c(t(genotypes[mk,])%*%betaG)

    tmp = (fk^a)*dnorm(fk/sigma)*(genotypes[k,gint]^c) - (fmk^a)*dnorm(fmk/sigma)*(genotypes[mk,gint]^c)

    tmp[(uindex[[u1]]==0)] = 0

    tmp2 = rep(0,length(y))
    for(u in 1:nu){
        cy = y[(uindex[[u]]==1)]
        cxe = as.matrix(xe[(uindex[[u]]==1),])
        cgenoprobs = genoprobs[[u]][,1]
        for (i in 1:dim(genotypes)[1]){
            fi = cy - alpha - cxe%*%betaE - c(t(genotypes[i,])%*%betaG)
            tmp2[(uindex[[u]]==1)] = tmp2[(uindex[[u]]==1)] + dnorm(fi/sigma)*cgenoprobs[i]
        }
    }
    return(tmp/tmp2)
}


chfun1dash = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,k,uindex,u1){
    nu = length(uindex)
    fk = y - alpha - xe%*%betaE - c(t(genotypes[k,])%*%betaG)
    mk = dim(genotypes)[1]
    fmk = y - alpha - xe%*%betaE - c(t(genotypes[mk,])%*%betaG)

    tmp = (fk^a)*dnorm(fk/sigma)*(genotypes[k,gint]^c)*genotypes[k,] - (fmk^a)*dnorm(fmk/sigma)*(genotypes[mk,gint]^c)*genotypes[mk,]

    tmp[(uindex[[u1]]==0)] = 0

    tmp2 = rep(0,length(y))
    for(u in 1:nu){
        cy = y[(uindex[[u]]==1)]
        cxe = as.matrix(xe[(uindex[[u]]==1),])
        cgenoprobs = genoprobs[[u]][,1]
        for (i in 1:dim(genotypes)[1]){
            fi = cy - alpha - cxe%*%betaE - c(t(genotypes[i,])%*%betaG)
            tmp2[(uindex[[u]]==1)] = tmp2[(uindex[[u]]==1)] + dnorm(fi/sigma)*cgenoprobs[i]
        }
    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}

chfun1 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,uindex){
    nu = length(uindex)
    tmp = matrix(0,ncol=dim(genotypes)[2],nrow = length(y))
    tmp2 = rep(0,length(y))
    for(u in 1:nu){
        cy = y[(uindex[[u]]==1)]
        cxe = as.matrix(xe[(uindex[[u]]==1),])
        cgenoprobs = genoprobs[[u]][,1]
        for (j in 1:dim(genotypes)[1]){
            fi = cy - alpha - cxe%*%betaE - c(t(genotypes[j,])%*%betaG)
            for(k in 1:dim(genotypes)[2]){
                tmp[(uindex[[u]]==1),k] = tmp[(uindex[[u]]==1),k] + ((fi^a)*dnorm(fi/sigma)*cgenoprobs[j]*(genotypes[j,gint]^c))*c(genotypes[j,k])
            }
            tmp2[(uindex[[u]]==1)] = tmp2[(uindex[[u]]==1)] + dnorm(fi/sigma)*cgenoprobs[j]
        }

    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}


chfun2 = function(y,xe,alpha,betaE,betaG,sigma,a,c,gint,genotypes,genoprobs,uindex){
    nu = length(uindex)
    tmp = matrix(0,ncol=dim(genotypes)[2]*dim(genotypes)[2],nrow = length(y))
    tmp2 = rep(0,length(y))
    for(u in 1:nu){
        cy = y[(uindex[[u]]==1)]
        cxe = as.matrix(xe[(uindex[[u]]==1),])
        cgenoprobs = genoprobs[[u]][,1]
        for (j in 1:dim(genotypes)[1]){
            fi = cy - alpha - cxe%*%betaE - c(t(genotypes[j,])%*%betaG)
            for(k in 1:dim(genotypes)[2]){
                tmp[(uindex[[u]]==1),k] = tmp[(uindex[[u]]==1),k] + ((fi^a)*dnorm(fi/sigma)*cgenoprobs[j]*(genotypes[j,gint]^c))*(as.vector(t(c(genotypes[j,]%*%t(genotypes[j,]))))[k])
            }
            tmp2[(uindex[[u]]==1)] = tmp2[(uindex[[u]]==1)] + dnorm(fi/sigma)*cgenoprobs[j]
        }

    }
    res = (apply(tmp, 2, "/", tmp2))
    return(res)
}
