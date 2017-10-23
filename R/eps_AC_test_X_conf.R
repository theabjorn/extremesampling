

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
    gm = apply(xg,2,meanimpute)

    ## v1 = (g_m^T(I-H)g_m) / sigma2  - For each individual
    # (Xe^T Xe)^-1
    xetxeinv = solve(crossprod(Xe,Xe))
    # G_m^T H G_m
    tmpMat = xetxeinv%*%crossprod(Xe,gm)
    v1 = (colSums(gm*gm) - colSums( crossprod(Xe,gm)*tmpMat ))/sigma2

    ## v2 =  sum_j n_j var(g|j) ( sigma2 - var(z_j)) / sigma4

    calcvar2 = function(g){
        tmp = 0
        for(u in 1:nu){
            uindex_all = (xec == ux[u,])
            gu = g[uindex_all]
            varg = var(gu,na.rm=TRUE)
            zu = z[uindex_all]
            n_ucc = sum(!is.na(gu))

            tmp = tmp + n_ucc*varg*(sigma2 - var(zu[!is.na(gu)]))
        }
        return(tmp/(sigma2^2))
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
