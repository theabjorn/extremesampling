



eps_CC_test_simple_EX = function(y,xe,xg,l,u){

    n = length(y)
    ng = dim(xg)[2]
    Xe = cbind(rep(1,n),xe)

    # Make one test statistics for ng genetic variants simultaneuosly

    fit = epsCC.loglikmax(cbind(y,xe),l,u)
    beta = fit[1:(length(fit)-1)]
    sigma = fit[length(fit)]
    sigma2 = sigma*sigma

    xbeta = Xe%*%beta
    f = y-xbeta
    zl = (l-xbeta)/sigma
    zu = (u-xbeta)/sigma

    h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
    h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
    h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
    h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

    a = c(1 - h1 - h0*h0)
    b = c(- h0 - h2 - h0*h1)
    c = 2 - c(h1 + h3 + h1*h1)

    I11_11 = crossprod(Xe,Xe*a) # t(Xe)%*%(diag(a)%*%Xe)

    I11_22 = sum(c)

    I11_12 = crossprod(Xe,b) # t(Xe)%*%(b)
    I11_21 = t(I11_12)

    I11 = rbind(cbind(I11_11, I11_12),
                cbind(I11_21, I11_22))

    I22 = crossprod(xg,xg*a)

    I12 = crossprod(xg,cbind(Xe*a,b))

    Sigma = (I22 - I12%*%tcrossprod(solve(I11),I12))/sigma2

    s = crossprod(xg,y-xbeta + sigma*h0)/sigma2

    t = crossprod(s,crossprod(ginv(Sigma),s))

    pvalue = pchisq(t,ng,lower.tail=FALSE)
    statistic = t
    result = list(statistic,pvalue)
    names(result) = c("statistic","p.value")
    return(result)
}
