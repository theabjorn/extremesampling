

###############################################################
# EPS ALL CASE
# Environmental covariates (xe) present in the null model
# No confounding
###############################################################

eps_AC_test_x = function(y,xe,xg){

    n = length(y)
    ng = dim(xg)[2]

    fit = lm(y~xe) # Fit under H0
    z = fit$residuals
    sigma2 = sum(z^2)/n
    sigma = sqrt(sigma2)
    Xe = cbind(rep(1,n),xe)

    meanimpute = function(g){
        g[is.na(g)] = mean(g,na.rm=TRUE)
        return(g)
    }
    gm = apply(xg,2,meanimpute)

    ## v1 = (g_m^T(I-H)g_m) / sigma2  - For each individual
    # (Xe^T Xe)^-1
    xetxeinv = solve(crossprod(Xe,Xe))
    # G_m^T H G_m
    tmpMat = xetxeinv%*%crossprod(Xe,gm)
    v1 = (colSums(gm*gm) - colSums( crossprod(Xe,gm)*tmpMat ))/sigma2

    ## v2 =  n var(g) (sigma2 - var(z)) / sigma4

    calcvar2 = function(g){
        varg = var(g,na.rm=TRUE)
        n_cc = sum(!is.na(g))
        return( n_cc*varg*(sigma2 - var(z[!is.na(g)]))/(sigma2^2) )
    }
    v2 = apply(xg,2,calcvar2)

    s = crossprod(gm,z)/sigma2
    Sigma = v1 - v2

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
