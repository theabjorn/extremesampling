

eps_AC_test_nox = function(y,xg){

    n = length(y)
    ng = dim(xg)[2]


    fit = lm(y~1) # Fit under H0
    z = fit$residuals
    sigma2 = sum(z^2)/n
    sigma = sqrt(sigma2)

    Xe = matrix(1,ncol=1,nrow = n)

    meanimpute = function(g){
        g[is.na(g)] = mean(g,na.rm=TRUE)
        return(g)
    }
    xgm = apply(xg,2,meanimpute)
    veg = crossprod(Xe,xgm)
    veeinvveg = solve(crossprod(Xe,Xe))%*%veg
    var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2
    calcvar2 = function(g){
        varg = var(g,na.rm=TRUE)
        n_ic = sum(is.na(g))
        n_cc = sum(!is.na(g))
        return(varg*(sum(z[is.na(g)]^2) - sigma2*n_ic + (1/n_cc)*(sum(z[is.na(g)]))^2)/(sigma2^2))
    }
    var2 = apply(xg,2,calcvar2)

    s = c(t(z)%*%xgm/sigma2)
    #s = c(t(z[!is.na(xg[,1])])%*%xg[!is.na(xg[,1]),]/sigma2)
    Sigma = var1 - var2
    t = (s*s)/Sigma
    pval = pchisq(t,1,lower.tail=FALSE)

    statistic = matrix(t,ncol = 1, nrow = ng)
    pvalue = matrix(pval,ncol = 1, nrow = ng)
    rownames(statistic) = colnames(xg)
    rownames(pvalue) =  colnames(xg)
    colnames(statistic) = "t"
    colnames(pvalue) = "p.value"
    result = list(statistic,pvalue)
    names(result) = c("statistic","p.value")
    return(result)
}


