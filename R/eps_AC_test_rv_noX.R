

###

eps_AC_test_rv_noX = function(y,xg,w){

    n = length(y)
    ng = dim(xg)[2]

    fit = lm(y~1) # Fit under H0
    z = fit$residuals
    sigma2 = sum(z^2)/n
    sigma = sqrt(sigma2)

    Xe = cbind(rep(1,n))

    meanimpute = function(g){
        g[is.na(g)] = mean(g,na.rm=TRUE)
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
    extreme = !is.na(xg[,1])
    n_cc = sum(extreme)
    v2 = n_cc*var(xgw,na.rm = TRUE)*(sigma2 - var(z[extreme]))/(sigma2^2)

    Sigma = v1 - v2

    s = crossprod(gmw,z)/sigma2

    return(list(s,Sigma))
}
