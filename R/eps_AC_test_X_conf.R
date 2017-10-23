

###############################################################
# EPS ALL CASE
# Environmental covariates (xe) present in the null model
# Confounding assumed
###############################################################

eps_AC_test_x_conf = function(y,xe,xec,xg){

    n = length(y)
    ng = dim(xg)[2]

    fit = lm(y~xe) # Fit under H0
    z = fit$residuals
    sigma2 = sum(z^2)/n
    sigma = sqrt(sigma2)

    Xe = cbind(rep(1,n),xe)

    ux = as.matrix(unique(xec))
    nu = dim(ux)[1]

    meanimpute = function(g){
        for(u in 1:nu){
            uindex_all = (xec == ux[u,])
            gu = g[uindex_all]
            gu[is.na(gu)] = mean(gu,na.rm=TRUE)
            g[uindex_all] = gu
        }
        return(g)
    }

    xgm = apply(xg,2,meanimpute)
    veg = crossprod(Xe,xgm)
    veeinvveg = solve(crossprod(Xe,Xe))%*%veg
    var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2

    calcvar2 = function(g){
        tmp = 0
        for(u in 1:nu){
            uindex_all = (xec == ux[u,])
            gu = g[uindex_all]
            varg = var(gu,na.rm=TRUE)
            zu = z[uindex_all]
            n_uic = sum(is.na(gu))
            n_ucc = sum(!is.na(gu))
            tmp = tmp + varg*(sum(zu[is.na(gu)]^2) - sigma2*n_uic + (1/n_ucc)*(sum(zu[is.na(gu)]))^2)/(sigma2^2)
        }
        return(tmp)
    }
    var2 = apply(xg,2,calcvar2)

    s = c(t(z)%*%xgm/sigma2)
    Sigma = var1 - var2
    t = (s*s)/Sigma
    pval = pchisq(t,1,lower.tail=FALSE)

    statistic = matrix(t,ncol = 1, nrow = ng)
    pvalue = matrix(pval,ncol = 1, nrow = ng)
    rownames(statistic) = colnames(xg)
    rownames(pvalue) = colnames(xg)
    colnames(statistic) = "t"
    colnames(pvalue) = "p.value"
    result = list(statistic,pvalue)
    names(result) = c("statistic","p.value")
    return(result)
}
