

###################################################################
# EPS-CC extreme-phenotype individuals only, covariates
###################################################################

eps_CC_test_ex = function(y,xe,xg,l,u){
    n = length(y)
    ng = dim(xg)[2]
    Xe = cbind(rep(1,n),xe)

    # Make ng test statistics, each based on one genetic covariate

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
    # I11 equal for all g_i (ne+1 x ne+1) (beta, sigma)

    colSums(a*xg*xg)

    I22 = colSums(a*xg*xg)   # t(g_i)%*%Diag(a)%*%g_i, one value for each g_i

    I21 = crossprod(xg,cbind(Xe*a,b))
    I12 = crossprod(cbind(Xe*a,b),xg)

    tmpMat = crossprod(I12,solve(I11))

    Sigma = (1/sigma2)*(I22 - rowSums(I21*tmpMat)) # one value for each g_i

    s = c(crossprod(xg,(y-xbeta+sigma*h0))/sigma2) # one value for each g_i

    t = s*s/Sigma  # one value for each g_i
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
