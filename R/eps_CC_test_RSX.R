
###################################################################
# EPS-CC extreme-phenotype individuals + random samples, covariates
###################################################################

eps_CC_test_rsx = function(y,xe,xg,l,u,randomindex){

    y_r = y[randomindex ==1]
    y_e = y[randomindex ==0]

    x_r = as.matrix(as.matrix(xe)[randomindex ==1,])
    x_e = as.matrix(as.matrix(xe)[randomindex ==0,])

    g_r = as.matrix(as.matrix(xg)[randomindex ==1,])
    g_e = as.matrix(as.matrix(xg)[randomindex ==0,])
    ng = dim(g_e)[2]

    nr = length(y_r)
    ne = length(y_e)

    fit = epsCC.loglikmax(cbind(y,xe),c(l,u),randomindex)
    alpha = fit[1]
    beta = fit[2:(length(fit)-1)]
    sigma = fit[length(fit)]
    sigma2 = sigma*sigma

    f_r = y_r-alpha - x_r%*%beta
    f_e = y_e-alpha - x_e%*%beta

    xbeta = x_e%*%beta
    zl = (l-alpha-xbeta)/sigma
    zu = (u-alpha-xbeta)/sigma

    h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
    h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
    h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
    h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

    a = c(1 - h1 - h0*h0)
    c = c(2*h1 - h3 - h1*h1)
    b = c(h0 - h2 - h0*h1)

    I11_11 = nr + sum(a)
    I11_22 = t(x_r)%*%x_r + t(x_e)%*%diag(a)%*%x_e

    I11_21 = colSums(x_r) + t(x_e)%*%a
    I11_12 = t(I11_21)

    I11_33 = sum(c - 3*c(h1)) + 2*ne + 2*nr
    I11_31 = sum(b-2*c(h0))
    I11_13 = t(I11_31)

    I11_23 = t(x_e)%*%(b-2*c(h0))
    I11_32 = t(I11_23)

    I11 = cbind(rbind(I11_11,I11_21,I11_31),
                rbind(I11_12,I11_22,I11_32),
                rbind(I11_13,I11_23,I11_33))

    statistic = matrix(NA,ncol = 1, nrow = ng)
    pvalue = matrix(NA,ncol = 1, nrow = ng)
    rownames(statistic) = colnames(xg)
    rownames(pvalue) = colnames(xg)
    colnames(statistic) = "t"
    colnames(pvalue) = "p.value"

    for(i in 1:ng){
        gei = g_e[,i]
        gri = g_r[,i]

        I22 = t(gri)%*%gri + t(gei)%*%(diag(a)%*%gei)

        I21_1 = sum(gri) + t(gei)%*%a

        I21_2 = t(gri)%*%x_r + t(gei)%*%diag(a)%*%x_e

        I21_3 = 2*sum(f_r*gri)/sigma + 2*sum(f_e*gei)/sigma + sum(b*gei)

        I21 = cbind(I21_1,I21_2,I21_3)
        I12 = t(I21)

        Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
        s = t(gri)%*%(y_r-alpha-x_r%*%beta)/sigma2 + t(gei)%*%(y_e-alpha-xbeta+sigma*h0)/sigma2

        t = s*s/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        statistic[i,] = c(t)
        pvalue[i,] = c(pval)
    }
    result = list(statistic,pvalue)
    names(result) = c("statistic","p.value")
    return(result)

}

