

eps_AC_test_rv_X_conf = function(y,xe,xec,xg,w){

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

    gmw = t(t(gm)*w)
    xgw = t(t(xg)*w)

    ## v1 = (G_m^T(I-H)G_m) / sigma2  - combined for matrix of genetic vars
    # (Xe^T Xe)^-1
    xetxeinv = solve(crossprod(Xe,Xe))
    # G_m^T H G_m
    gmthgm = crossprod(gmw,Xe)%*%xetxeinv%*%crossprod(Xe,gmw)
    v1 = (crossprod(gmw,gmw) - gmthgm )/sigma2

    ## v2 =  n var(g) (sigma2 - var(z)) / sigma4

    tmp = 0*1
    for(u in 1:nu){
        uindex_all = (xec == ux[u,])
        gu = xgw[uindex_all,]
        varg = var(gu,na.rm=TRUE)
        zu = z[uindex_all]
        n_ucc = sum(!is.na(gu[,1]))
        tmp = tmp + n_ucc*varg*(sigma2 - var(zu[!is.na(gu[,1])]))
    }
    v2 = tmp/(sigma2^2)

    Sigma = v1 - v2

    s = crossprod(gmw,z)/sigma2

    return(list(s,Sigma))
}
