epsCC.rv.test.naive = function(epsdata0,covariates0,RV,isx,l,u,rsample, randomindex){

    if(!rsample){
        ###################################################################
        # No random sample, extreme-phenotype individuals only
        ###################################################################
        message("EPS complete-case analysis with no random samples")

        y = epsdata0[,1]
        n = length(y)

        if(sum((y>l & y<u))>0){stop("Incorrect data format")}

        g = as.matrix(RV)
        ng = dim(g)[2]

        modeldata = cbind(epsdata0[,1],covariates0)

        if(isx){
            ###############################################################
            # Covariates present in the null model
            ###############################################################
            x = as.matrix(covariates0)

            fit = epsCC.loglikmax(modeldata,c(l,u)) # Fit under H0
            alpha = fit[1]
            beta = fit[2:(length(fit)-1)]
            sigma = fit[length(fit)]
            sigma2 = sigma*sigma

            f = y-alpha - x%*%beta

            xbeta = x%*%beta
            zl = (l-alpha-xbeta)/sigma
            zu = (u-alpha-xbeta)/sigma

            h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
            h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
            h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
            h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

            a = c(1 - h1 - h0*h0)
            c = c(2*h1 - h3 - h1*h1)
            b = c(h0 - h2 - h0*h1)

            I11_11 = sum(a)
            I11_22 = t(x)%*%(diag(a)%*%x)
            I11_33 = sum(c - 3*c(h1)) + 2*n

            I11_21 = t(x)%*%a; I11_12 = t(I11_21)
            I11_31 = sum(b-2*c(h0)); I11_13 = t(I11_31)
            I11_23 = t(x)%*%(b-2*c(h0)); I11_32 = t(I11_23)

            I11 = cbind(rbind(I11_11,I11_21,I11_31),
                        rbind(I11_12,I11_22,I11_32),
                        rbind(I11_13,I11_23,I11_33))

            I22 = t(g)%*%(diag(a)%*%g)

            I21_1 = t(g)%*%a
            I21_2 = t(g)%*%(diag(a)%*%x)
            I21_3 = 2*t(g)%*%f/sigma + t(g)%*%b

            I21 = cbind(I21_1,I21_2,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y-alpha - xbeta +sigma*h0)%*%g/sigma2
            t = s%*%ginv(Sigma)%*%t(s)
            pvalue = pchisq(t,ng,lower.tail=FALSE)
            statistic = t
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }else{
            ###############################################################
            # No covariates present in the null model
            ###############################################################
            fit = epsCC.loglikmax(modeldata,c(l,u)) # Fit under H0
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
            c = c(2*h1 - h3 - h1*h1)
            b = c(h0 - h2 - h0*h1)

            I11_11 = n*a
            I11_33 = n*(c - 3*c(h1)) + 2*n

            I11_31 = n*(b-2*c(h0)); I11_13 = t(I11_31)

            I11 = cbind(rbind(I11_11,I11_31),
                        rbind(I11_13,I11_33))

            I22 = a*t(g)%*%g

            I21_1 = a*colSums(g)
            I21_3 = 2*t(g)%*%f/sigma + b*colSums(g)

            I21 = cbind(I21_1,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y-alpha+sigma*h0)%*%g/sigma2
            t = s%*%ginv(Sigma)%*%t(s)
            pvalue = pchisq(t,ng,lower.tail=FALSE)
            statistic = t
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }
    }else{
        ###################################################################
        # Extremes + random sample
        ###################################################################
        message("EPS complete-case analysis with random samples")

        y_r = epsdata0[,1][randomindex ==1]
        y_e = epsdata0[,1][randomindex ==0]
        nr = length(y_r)
        ne = length(y_e)

        gr = as.matrix(as.matrix(RV)[randomindex ==1,])
        ge = as.matrix(as.matrix(RV)[randomindex ==0,])
        ng = dim(ge)[2]

        modeldata = epsdata0

        if(isx){
            ###############################################################
            # Covariates present in the null model
            ###############################################################
            x_r = as.matrix(as.matrix(covariates0)[randomindex ==1,])
            x_e = as.matrix(as.matrix(covariates0)[randomindex ==0,])

            fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
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

            I22 = t(gr)%*%gr + t(ge)%*%(diag(a)%*%ge)

            I21_1 = colSums(gr) + t(ge)%*%a
            I21_2 = t(gr)%*%x_r + t(ge)%*%diag(a)%*%x_e

            I21_3 = 2*t(gr)%*%f_r/sigma + 2*t(ge)%*%f_e/sigma + t(ge)%*%b

            I21 = cbind(I21_1,I21_2,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y_r-alpha - x_r%*%beta)%*%gr/sigma2 + t(y_e-alpha - xbeta +sigma*h0)%*%ge/sigma2

            t = s%*%ginv(Sigma)%*%t(s)
            pvalue = pchisq(t,ng,lower.tail=FALSE)
            statistic = t
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }else{
            ##############################################################
            # No covariates present in the null model
            ##############################################################
            fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
            alpha = fit[1]
            sigma = fit[length(fit)]
            sigma2 = sigma*sigma

            f_r = y_r-alpha
            f_e = y_e-alpha

            zl = (l-alpha)/sigma
            zu = (u-alpha)/sigma

            h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
            h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
            h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
            h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

            a = c(1 - h1 - h0*h0)
            c = c(2*h1 - h3 - h1*h1)
            b = c(h0 - h2 - h0*h1)

            I11_11 = nr + a*ne

            I11_33 = ne*(c - 3*c(h1)) + 2*ne + 2*nr
            I11_31 = ne*(b-2*c(h0))
            I11_13 = t(I11_31)

            I11 = cbind(rbind(I11_11,I11_31),
                        rbind(I11_13,I11_33))

            I22 = t(gr)%*%gr + a*t(ge)%*%ge

            I21_1 = colSums(gr) + a*colSums(ge)
            I21_3 = 2*t(gr)%*%f_r/sigma + 2*t(ge)%*%f_e/sigma + b*colSums(ge)

            I21 = cbind(I21_1,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y_r-alpha)%*%gr/sigma2 + t(y_e-alpha+sigma*h0)%*%ge/sigma2

            t = s%*%ginv(Sigma)%*%t(s)
            pvalue = pchisq(t,ng,lower.tail=FALSE)
            statistic = t
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }
    }
}



epsCC.rv.test.lmm = function(epsdata0,covariates0,RV,isx,l,u,rsample,randomindex,weights){

    Bmat = diag(weights)

    if(!rsample){
        ###################################################################
        # No random sample, extreme-phenotype individuals only
        ###################################################################
        message("EPS complete-case analysis with no random samples")

        y = epsdata0[,1]
        n = length(y)

        if(sum((y>l & y<u))>0){stop("Incorrect data format")}

        g = as.matrix(RV)
        ng = dim(g)[2]

        modeldata = cbind(epsdata0[,1],covariates0)

        if(isx){
            ###############################################################
            # Covariates present in the null model
            ###############################################################
            x = as.matrix(covariates0)

            fit = epsCC.loglikmax(modeldata,c(l,u)) # Fit under H0
            alpha = fit[1]
            beta = fit[2:(length(fit)-1)]
            sigma = fit[length(fit)]
            sigma2 = sigma*sigma

            f = y-alpha - x%*%beta

            xbeta = x%*%beta
            zl = (l-alpha-xbeta)/sigma
            zu = (u-alpha-xbeta)/sigma

            h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
            h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
            h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
            h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/(1-pnorm(zu)+pnorm(zl))

            a = c(1 - h1 - h0*h0)
            c = c(2*h1 - h3 - h1*h1)
            b = c(h0 - h2 - h0*h1)

            I11_11 = sum(a)
            I11_22 = t(x)%*%(diag(a)%*%x)
            I11_33 = sum(c - 3*c(h1)) + 2*n

            I11_21 = t(x)%*%a; I11_12 = t(I11_21)
            I11_31 = sum(b-2*c(h0)); I11_13 = t(I11_31)
            I11_23 = t(x)%*%(b-2*c(h0)); I11_32 = t(I11_23)

            I11 = cbind(rbind(I11_11,I11_21,I11_31),
                        rbind(I11_12,I11_22,I11_32),
                        rbind(I11_13,I11_23,I11_33))

            I22 = t(g)%*%(diag(a)%*%g)

            I21_1 = t(g)%*%a
            I21_2 = t(g)%*%(diag(a)%*%x)
            I21_3 = 2*t(g)%*%f/sigma + t(g)%*%b

            I21 = cbind(I21_1,I21_2,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y-alpha - xbeta +sigma*h0)%*%g/sigma2

            sqrtSigma = sqrtm(Sigma)
            eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
            d = eigen(eigenmat)$values

            qdist = rep(0,1000000)
            for(c in 1:length(d)){
                qdist = qdist + d[c]*rchisq(1000000,1)
            }

            #critval = quantile(qdist,probs = 0.95)

            teststat = c(s%*%Bmat%*%t(s))

            pvalue = sum(qdist > teststat)/1000000

            result = list(teststat,pvalue,d)
            names(result) = c("statistic","p.value","eigen")
            return(result)
        }else{
            ###############################################################
            # No covariates present in the null model - use expected info matrix
            ###############################################################
            fit = epsCC.loglikmax(modeldata,c(l,u)) # Fit under H0
            alpha = fit[1]
            sigma = fit[length(fit)]
            sigma2 = sigma*sigma

            zl = (l-alpha)/sigma
            zu = (u-alpha)/sigma

            h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
            h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))

            a = c(1 - h1 - h0*h0)

            O = matrix(1,nrow = n,ncol = 1); H = O%*%ginv(t(O)%*%O)%*%t(O); I = diag(1,n)

            Sigma = (a/sigma2)*t(g)%*%(I-H)%*%g
            s = t(y-mean(y))%*%g/sigma2

            sqrtSigma = sqrtm(Sigma)
            eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
            d = eigen(eigenmat)$values

            qdist = rep(0,1000000)
            for(c in 1:length(d)){
                qdist = qdist + d[c]*rchisq(1000000,1)
            }

            #critval = quantile(qdist,probs = 0.95)

            teststat = c(s%*%Bmat%*%t(s))

            pvalue = sum(qdist > teststat)/1000000

            result = list(teststat,pvalue,d)
            names(result) = c("statistic","p.value","eigen")
            return(result)
        }
    }else{
        ###################################################################
        # Extremes + random sample
        ###################################################################
        message("EPS complete-case analysis with random samples")

        y_r = epsdata0[,1][randomindex ==1]
        y_e = epsdata0[,1][randomindex ==0]
        nr = length(y_r)
        ne = length(y_e)
        n = nr + ne

        g = as.matrix(RV)
        gr = g[randomindex ==1,]
        ge = g[randomindex ==0,]
        ng = dim(g)[2]

        modeldata = epsdata0

        if(isx){
            ###############################################################
            # Covariates present in the null model
            ###############################################################
            x_r = as.matrix(as.matrix(covariates0)[randomindex ==1,])
            x_e = as.matrix(as.matrix(covariates0)[randomindex ==0,])

            fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
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

            I22 = t(gr)%*%gr + t(ge)%*%(diag(a)%*%ge)

            I21_1 = colSums(gr) + t(ge)%*%a
            I21_2 = t(gr)%*%x_r + t(ge)%*%diag(a)%*%x_e

            I21_3 = 2*t(gr)%*%f_r/sigma + 2*t(ge)%*%f_e/sigma + t(ge)%*%b

            I21 = cbind(I21_1,I21_2,I21_3); I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = t(y_r-alpha - x_r%*%beta)%*%gr/sigma2 + t(y_e-alpha - xbeta +sigma*h0)%*%ge/sigma2

            sqrtSigma = sqrtm(Sigma)
            eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
            d = eigen(eigenmat)$values

            qdist = rep(0,1000000)
            for(c in 1:length(d)){
                qdist = qdist + d[c]*rchisq(1000000,1)
            }

            #critval = quantile(qdist,probs = 0.95)

            teststat = c(s%*%Bmat%*%t(s))

            pvalue = sum(qdist > teststat)/10000000

            result = list(teststat,pvalue,d)
            names(result) = c("statistic","p.value","eigen")
            return(result)
        }else{
            ##############################################################
            # No covariates present in the null model - use expected info
            ##############################################################
            fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
            alpha = fit[1]
            sigma = fit[length(fit)]
            sigma2 = sigma*sigma

            f_r = y_r-alpha
            f_e = y_e-alpha

            zl = (l-alpha)/sigma
            zu = (u-alpha)/sigma

            h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
            h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))

            a2 = c(1- (ne/n)*h1 - (ne/n)*h0*h0)

            O = matrix(1,nrow = n,ncol = 1); H = O%*%ginv(t(O)%*%O)%*%t(O); I = diag(1,n)

            Sigma = (a2/sigma2)*t(g)%*%(I-H)%*%g
            s = t(y-mean(y))%*%g/sigma2

            sqrtSigma = sqrtm(Sigma)
            eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
            d = eigen(eigenmat)$values

            qdist = rep(0,1000000)
            for(c in 1:length(d)){
                qdist = qdist + d[c]*rchisq(1000000,1)
            }

            #critval = quantile(qdist,probs = 0.95)

            teststat = c(s%*%Bmat%*%t(s))

            pvalue = sum(qdist > teststat)/10000000

            result = list(teststat,pvalue,d)
            names(result) = c("statistic","p.value","eigen")
            return(result)
        }
    }
}
