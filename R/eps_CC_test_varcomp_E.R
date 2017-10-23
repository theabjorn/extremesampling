

eps_CC_test_varcomp_E = function(y,xg,l,u,weights){

    n = length(y)
    ng = dim(xg)[2]

    fit = epsCC.loglikmax(cbind(y),l,u) # Fit under H0
    alpha = fit[1]
    sigma = fit[length(fit)]
    sigma2 = sigma*sigma

    f = y-alpha
    zl = (l-alpha)/sigma
    zu = (u-alpha)/sigma

    h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
    h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
    h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
    h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

    a = c(1 - h1 - h0*h0)
    b = c(- h0 - h2 - h0*h1)
    c = 2 - c(h1 + h3 + h1*h1)

    I11_11 = n*a
    I11_33 = n*c

    I11_31 = n*b; I11_13 = t(I11_31)

    I11 = cbind(rbind(I11_11,I11_31),
                rbind(I11_13,I11_33))

    I22 = a*crossprod(xg,xg)

    I12 = crossprod(xg,cbind(rep(a,n),rep(b,n)))

    Sigma = (I22 - I12%*%tcrossprod(solve(I11),I12))/sigma2

    s = crossprod(xg,y-alpha + sigma*h0)/sigma2

    d = eigen(Sigma)$values

    teststat = crossprod(s)
    #pvalue = sum(qdist > teststat)/1000000
    pvalue = davies(teststat,d,acc = 0.00005)$Qq

    result = list(teststat,pvalue,d)
    names(result) = c("statistic","p.value","eigen")
    return(result)
}
