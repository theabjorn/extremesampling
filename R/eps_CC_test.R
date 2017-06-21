#' @title Test for associations under the EPS complete-case design
#' @description
#' \code{epsCC.test} performs a score test for common genetic variants
#' under the EPS complete-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param SNP a matrix of genetic variants to be tested against the null
#' @param cutoffs a vector \code{c(l,u)} of the lower and upper cut-offs used
#' for extreme sampling
#' @param randomindex a binary vector that indicates if samples are random
#' or extreme
#' @return \code{epsCC.test} returns for each genetic variant:
#' \item{statistic}{the score test statistic}
#' \item{p.value}{the P-value}
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe, which describes a regression model, y=a+be*xe+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The variables are taken from the environment that the
#' function is called from.
#' The null hypothesis bg=0 is tested for the model y=a+be*xe+bg*xg+e.
#' The covariate \code{xg} is a SNP (single-nucleotide polymorphism).
#' Both xe and xg can be matrices.
#'
#' The EPS complete-case design is such that data is only available
#' for individuals with high and low values of the phenotype \code{y}, and
#' potentialy some randomly sampled individuals. The cut-offs \code{l} and
#' \code{u} that specify the sampling must be specified
#' in the \code{cutoffs} argument.
#'
#' @import MASS stats
#' @export
#' @examples
#' N = 5000 # Number of individuals in a population
#' xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
#' xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
#' xg1 = sample(c(0,1,2),N,c(0.4,0.3,0.3), replace = TRUE) # SNP
#' xg2 = sample(c(0,1,2),N,c(0.5,0.3,0.2), replace = TRUE) # SNP
#' # Model parameters
#' a = 50; be1 = 5; be2 = 8; bg1 = 0.3; bg2 = 0.6; sigma = 2
#' # Generate response y
#' y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 + bg2*xg2, sd = sigma)
#' # Identify extremes, here upper and lower 25% of population
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' # Create the EPS-only data set
#' y = y[extreme]
#' xe1 = xe1[extreme]
#' xe2 = xe2[extreme]
#' xg1 = xg1[extreme]
#' xg2 = xg2[extreme]
#' xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#'
#' epsCC.test(y~xe,SNP=xg,cutoffs = c(l,u))
#' epsCC.test(y~xe,SNP=xg,cutoffs = c(l,u),onebyone = FALSE)

epsCC.test = function(nullmodel,SNP,cutoffs,randomindex){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(length(cutoffs) != 2){stop("Invalid cutoffs vector given")}

    rsample = TRUE
    if(missing(randomindex)){
        rsample = FALSE
    }else if(sum(randomindex)==0){
        rsample = FALSE
    }

    if(is.null(colnames(SNP))){
        SNP = as.matrix(SNP)
        colnames(SNP) = paste0("SNP",1:dim(SNP)[2])}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    isx = (dim(covariates0)[2]>0) # there are covariates in the null model

    if(sum(is.na(epsdata0) > 1)){
        stop("NA values in the data not allowed")}

    l = min(cutoffs)
    u = max(cutoffs)

    if(!rsample){
        ###################################################################
        # No random sample, extreme-phenotype individuals only
        ###################################################################
        message("EPS complete-case analysis with no random samples")

        y = epsdata0[,1]
        n = length(y)

        if(sum((y>l & y<u))>0){
            stop("Incorrect data format")}

        totest = colnames(SNP)
        g = as.matrix(SNP)
        ng = dim(g)[2]

        modeldata = cbind(epsdata0[,1],covariates0)

        ###################################################################
        # Test genetic covariates (one at a time)
        ###################################################################
        statistic = matrix(NA,ncol = 1, nrow = ng)
        pvalue = matrix(NA,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"

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

            xbeta = x%*%beta
            f = y-alpha-xbeta
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

            I11_21 = t(x)%*%a
            I11_12 = t(I11_21)

            I11_33 = sum(c - 3*c(h1)) + 2*n

            I11_31 = sum(b-2*c(h0))
            I11_13 = t(I11_31)

            I11_23 = t(x)%*%(b-2*c(h0))
            I11_32 = t(I11_23)

            I11 = cbind(rbind(I11_11,I11_21,I11_31),
                        rbind(I11_12,I11_22,I11_32),
                        rbind(I11_13,I11_23,I11_33))

            for(i in 1:ng){
                gi = g[,i]

                I22 = t(gi)%*%(diag(a)%*%gi)

                I21_1 = t(gi)%*%a
                I21_2 = t(gi)%*%diag(a)%*%x

                I21_3 = 2*sum(f*gi)/sigma + sum(b*gi)

                I21 = cbind(I21_1,I21_2,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = t(gi)%*%(y-alpha-xbeta +sigma*h0)/sigma2
                t = s*s/Sigma
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[i,] = c(t)
                pvalue[i,] = c(pval)
            }
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }else{
            ##############################################################
            # No covariates present in the null model
            ##############################################################
            # = score test for random sample
            fit = epsCC.loglikmax(as.matrix(y),c(l,u)) # Fit under H0
            alpha = fit[1]
            sigma = fit[2]
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

            I11_31 = n*(b-2*c(h0))
            I11_13 = t(I11_31)

            I11 = cbind(rbind(I11_11,I11_31),
                        rbind(I11_13,I11_33))

            for(i in 1:ng){
                gi = g[,i]

                I22 = a*t(gi)%*%gi

                I21_1 = a*sum(gi)
                I21_3 = 2*sum(f*gi)/sigma + b*sum(gi)

                I21 = cbind(I21_1,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = t(gi)%*%(y-alpha+sigma*h0)/sigma2
                t = s*s/Sigma
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[i,] = c(t)
                pvalue[i,] = c(pval)
            }
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

        totest = colnames(SNP)
        g_r = as.matrix(as.matrix(SNP)[randomindex ==1,])
        g_e = as.matrix(as.matrix(SNP)[randomindex ==0,])
        ng = dim(g_e)[2]

        modeldata = epsdata0

        ###################################################################
        # Test genetic covariates (one at a time)
        ###################################################################
        statistic = matrix(NA,ncol = 1, nrow = ng)
        pvalue = matrix(NA,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"

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
        }else{
            ##############################################################
            # No covariates present in the null model
            ##############################################################
            fit = epsonlyloglikmax(as.matrix(y),c(l,u),randomindex) # Fit under H0
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

            I11_11 = nr + ne*a

            I11_33 = ne*(c - 3*c(h1)) + 2*ne + 2*nr
            I11_31 = ne*(b-2*c(h0))
            I11_13 = t(I11_31)

            I11 = cbind(rbind(I11_11,I11_31),
                        rbind(I11_13,I11_33))

            for(i in 1:ng){
                gei = g_e[,i]
                gri = g_r[,i]

                I22 = t(gri)%*%gri + a*t(gei)%*%gei

                I21_1 = sum(gri) + a*sum(gei)

                I21_3 = 2*sum(f_r*gri)/sigma + 2*sum(f_e*gei)/sigma + sum(b*gei)

                I21 = cbind(I21_1,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = t(gri)%*%(y_r-alpha)/sigma2 + t(gei)%*%(y_e-alpha+sigma*h0)/sigma2

                t = s*s/Sigma
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[i,] = c(t)
                pvalue[i,] = c(pval)
            }
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }
    }
}
