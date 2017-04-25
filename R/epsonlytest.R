#' @title Test for associations under the EPS-only design
#' @description
#' \code{epsonly.test} performs a score test for genetic variables
#' in the EPS-only design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param SNP a matrix of variables to be tested against the null
#' @param cutoffs a vector \code{c(l,u)} of the lower and upper cut-offs used
#' for extreme sampling
#' @param randomindex a vector which indicates if observations are from a random
#' or extreme sample
#' @param onebyone \code{TRUE} if genetic variables should be tested one by one
#' for inclusion in the model, default \code{TRUE}
#' @return \code{epsonly.test} returns
#' \item{statistic}{the value of the score test statistic}
#' \item{parameter}{the degrees of freedom of the statistic}
#' \item{p.value}{the P-value for the test}
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe, which describes a regression model, y=a+be*xe+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The variables are taken from the environment that the
#' function is called from.
#' The test considers a regression model y=a+be*xe+bg*xg+e, where xg is a
#' genetic covariate, and where bg=0 under the null hypothesis. The output
#' of the function gives the test statistic and p-value for the test of
#' H0: bg=0.
#' The covariate \code{xg} is a SNP (single-nucleotide polymorphism).
#' Both xe and xg can be matrices.
#'
#' The EPS-only design is such that data is only available
#' for individuals with high and low values of the phenotype \code{y}. The
#' cut-offs \code{l} and \code{u} that specify the sampling must be specified
#' in the \code{cutoffs} argument.
#'
#' Thus, the data set must consist of observations only for extreme individuals
#' \code{y > u} or \code{y < l} where \code{y} is the response (phenotype)
#' of the linear regression model.
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
#' epsonly.test(y~xe,SNP=xg,cutoffs = c(l,u))
#' epsonly.test(y~xe,SNP=xg,cutoffs = c(l,u),onebyone = FALSE)

epsonly.test = function(nullmodel,SNP,cutoffs,randomindex,onebyone = TRUE){
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

        if(sum(is.na(epsdata0) > 1)){
            stop("NA values in the data not allowed")}

        totest = colnames(SNP)
        g = as.matrix(SNP)

        modeldata = cbind(epsdata0[,1],covariates0)

        isx = (dim(covariates0)[2]>0) # there are covariates in the null model

        ng = dim(g)[2]
        n = dim(epsdata0)[1]
        y = epsdata0[,1]

        if(onebyone){
            ###################################################################
            # Test genetic covariates one by one (GWAS)
            ###################################################################
            statistic = matrix(NA,ncol = 1, nrow = ng)
            parameter = matrix(NA,ncol = 1, nrow = ng)
            pvalue = matrix(NA,ncol = 1, nrow = ng)
            rownames(statistic) = totest
            rownames(parameter) = totest
            rownames(pvalue) = totest
            colnames(statistic) = "t"
            colnames(parameter) = "d.o.f"
            colnames(pvalue) = "p.value"

            if(isx){
                ###############################################################
                # Covariates present in the null model
                ###############################################################
                x = as.matrix(covariates0)

                fit = epsonlyloglikmax(modeldata,c(l,u)) # Fit under H0
                alpha = fit[1]
                beta = fit[2:(length(fit)-1)]
                sigma = fit[length(fit)]
                sigma2 = sigma*sigma
                xbeta = x%*%beta
                zl = (l-alpha-xbeta)/sigma
                zu = (u-alpha-xbeta)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = sum(a)
                I11_22 = t(x)%*%(diag(a[,1])%*%x)

                I11_21 = t(t(a)%*%x)
                I11_12 = t(a)%*%x

                I11_33 = sum(c)
                I11_31 = sum(b)
                I11_13 = sum(b)
                I11_23 = t(t(b)%*%x)
                I11_32 = t(b)%*%x

                I11 = cbind(rbind(I11_11,I11_21,I11_31),
                            rbind(I11_12,I11_22,I11_32),
                            rbind(I11_13,I11_23,I11_33))

                for(i in 1:ng){
                    gi = g[,i]
                    I22 = sum(gi*gi*a[,1])

                    I21_1 = sum(a[,1]*gi)
                    I21_2 = t(t(x)%*%(a[,1]*gi))
                    I21_3 = sum(b[,1]*gi)
                    I21 = cbind(I21_1,I21_2,I21_3)
                    I12 = t(I21)

                    Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                    s = sum((y-alpha - xbeta +sigma*h0)*gi)/sigma2
                    t = s*s/Sigma
                    pval = pchisq(t,1,lower.tail=FALSE)

                    statistic[i,] = t
                    parameter[i,] = 1
                    pvalue[i,] = pval
                }
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }else{
                ##############################################################
                # No covariates present in the null model
                ##############################################################

                # score test for random sample

                y = epsdata0[,1]

                alpha = mean(y) # fit under H0
                sigma = sd(y)
                sigma2 = sigma*sigma

                I11_11 = n
                I11_33 = 2*n

                I11_31 = 0
                I11_13 = 0
                I11 = cbind(rbind(I11_11,I11_31),
                            rbind(I11_13,I11_33))

                for(i in 1:ng){
                    gi = g[,i]
                    I22 = t(gi)%*%gi
                    I21_1 = sum(gi)
                    I21_3 = 0
                    I21 = cbind(I21_1,I21_3)
                    Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                    s = (1/sigma2)*sum((y-alpha)*gi)
                    t = s*s/Sigma
                    pval = pchisq(t,1,lower.tail=FALSE)

                    statistic[i,] = t
                    parameter[i,] = 1
                    pvalue[i,] = pval
                }
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
        }else{
            ###################################################################
            # Test all genetic covariates at a time (canditate SNP)
            ###################################################################
            if(ng > 10){
                stop("Do not test more than 10 SNPs simultaneously,
                     choose onebyone = TRUE")
            }
            statistic = matrix(NA,ncol = 1, nrow = 1)
            parameter = matrix(NA,ncol = 1, nrow = 1)
            pvalue = matrix(NA,ncol = 1, nrow = 1)
            colnames(statistic) = "t"
            colnames(parameter) = "d.o.f"
            colnames(pvalue) = "p.value"
            if(isx){
                ###############################################################
                # Covariates present in the null model
                ###############################################################
                x = covariates0

                fit = epsonlyloglikmax(modeldata,c(l,u))
                alpha = fit[1]
                beta = fit[2:(length(fit)-1)]
                sigma = fit[length(fit)]
                sigma2 = sigma*sigma

                xbeta = x%*%beta
                zl = (l-alpha-xbeta)/sigma
                zu = (u-alpha-xbeta)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = sum(a)
                I11_22 = t(x)%*%(diag(a[,1])%*%x)

                I11_21 = t(t(a)%*%x)
                I11_12 = t(a)%*%x

                I11_33 = sum(c)
                I11_31 = sum(b)
                I11_13 = sum(b)
                I11_23 = t(t(b)%*%x)
                I11_32 = t(b)%*%x

                I11 = cbind(rbind(I11_11,I11_21,I11_31),
                            rbind(I11_12,I11_22,I11_32),
                            rbind(I11_13,I11_23,I11_33))

                I22 = t(g)%*%(diag(a[,1])%*%g)

                I21_1 = t(t(a)%*%g)
                I21_2 = t(t(x)%*%(diag(a[,1])%*%g))
                I21_3 = t(t(b)%*%g)
                I21 = cbind(I21_1,I21_2,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = t(y-alpha - xbeta +sigma*h0)%*%g/sigma2
                t = s%*%ginv(Sigma)%*%t(s)
                pvalue[1,1] = pchisq(t,ng,lower.tail=FALSE)
                statistic[1,1] = t
                parameter[1,1] = ng
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }else{
                ###############################################################
                # No covariates present in the null model
                ###############################################################
                y = epsdata0[,1]

                alpha = mean(y)
                sigma = sd(y)
                sigma2 = sigma*sigma

                I11_11 = n
                I11_33 = 2*n
                I11_31 = 0
                I11_13 = 0
                I11 = cbind(rbind(I11_11,I11_31),
                            rbind(I11_13,I11_33))
                I22 = t(g)%*%(g)

                I21_1 = matrix(colSums(g))
                I21_3 = matrix(0,nrow=ng,ncol=1)
                I21 = cbind(I21_1,I21_3)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = t(y-alpha)%*%g/sigma2
                t = s%*%ginv(Sigma)%*%t(s)
                pvalue[1,1] = pchisq(t,ng,lower.tail=FALSE)
                statistic[1,1] = t
                parameter[1,1] = ng
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
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

        isx = (dim(covariates0)[2]>0) # there are covariates in the null model

        if(onebyone){
            ###################################################################
            # Test genetic covariates one by one (GWAS)
            ###################################################################
            statistic = matrix(NA,ncol = 1, nrow = ng)
            parameter = matrix(NA,ncol = 1, nrow = ng)
            pvalue = matrix(NA,ncol = 1, nrow = ng)
            rownames(statistic) = totest
            rownames(parameter) = totest
            rownames(pvalue) = totest
            colnames(statistic) = "t"
            colnames(parameter) = "d.o.f"
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

                xbeta = x_e%*%beta
                zl = (l-alpha-xbeta)/sigma
                zu = (u-alpha-xbeta)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = nr + sum(a)
                I11_22 = t(x_r)%*%x_r + t(x_e)%*%(diag(a[,1])%*%x_e)

                I11_21 = t(t(colSums(x_r))) + t(t(a)%*%x_e)
                I11_12 = t(I11_21)

                I11_33 = 2*nr + sum(c)
                I11_31 = sum(b)
                I11_13 = sum(b)
                I11_23 = t(t(b)%*%x_e)
                I11_32 = t(b)%*%x_e

                I11 = cbind(rbind(I11_11,I11_21,I11_31),
                            rbind(I11_12,I11_22,I11_32),
                            rbind(I11_13,I11_23,I11_33))

                for(i in 1:ng){
                    gi_r = g_r[,i]
                    gi_e = g_e[,i]

                    I22 = sum(gi_r*gi_r) + sum(gi_e*gi_e*a[,1])

                    I21_1 = sum(gi_r) + sum(a[,1]*gi_e)
                    I21_2 = t(t(x_r)%*%(gi_r)) + t(t(x_e)%*%(a[,1]*gi_e))
                    I21_3 = sum(gi_r) + sum(b[,1]*gi_e)
                    I21 = cbind(I21_1,I21_2,I21_3)
                    I12 = t(I21)

                    Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                    s = (sum((y_r - alpha - x_r%*%beta)*gi_r) + sum((y_e-alpha - xbeta +sigma*h0)*gi_e))/sigma2
                    t = s*s/Sigma
                    pval = pchisq(t,1,lower.tail=FALSE)

                    statistic[i,] = t
                    parameter[i,] = 1
                    pvalue[i,] = pval
                }
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }else{
                ##############################################################
                # No covariates present in the null model
                ##############################################################
                fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
                alpha = fit[1]
                sigma = fit[length(fit)]
                sigma2 = sigma*sigma

                zl = (l-alpha)/sigma
                zu = (u-alpha)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = nr + ne*a

                I11_33 = 2*nr + ne*c
                I11_31 = ne*b
                I11_13 = ne*b

                I11 = cbind(rbind(I11_11,I11_31),
                            rbind(I11_13,I11_33))

                for(i in 1:ng){
                    gi_r = g_r[,i]
                    gi_e = g_e[,i]

                    I22 = sum(gi_r*gi_r) + sum(gi_e*gi_e*a)

                    I21_1 = sum(gi_r) + sum(a*gi_e)
                    I21_3 = sum(gi_r) + sum(b*gi_e)
                    I21 = cbind(I21_1,I21_3)
                    I12 = t(I21)

                    Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                    s = (sum((y_r - alpha)*gi_r) + sum((y_e-alpha +sigma*h0)*gi_e))/sigma2
                    t = s*s/Sigma
                    pval = pchisq(t,1,lower.tail=FALSE)

                    statistic[i,] = t
                    parameter[i,] = 1
                    pvalue[i,] = pval
                }
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
        }else{
            ###################################################################
            # Test all genetic covariates at a time (canditate SNP)
            ###################################################################
            if(ng > 10){
                stop("Do not test more than 10 SNPs simultaneously,
                     choose onebyone = TRUE")
            }
            statistic = matrix(NA,ncol = 1, nrow = 1)
            parameter = matrix(NA,ncol = 1, nrow = 1)
            pvalue = matrix(NA,ncol = 1, nrow = 1)
            colnames(statistic) = "t"
            colnames(parameter) = "d.o.f"
            colnames(pvalue) = "p.value"
            if(isx){
                ###############################################################
                # Covariates present in the null model
                ###############################################################
                x_r = as.matrix(covariates0)[randomindex ==1,]
                x_e = as.matrix(covariates0)[randomindex ==0,]

                fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
                alpha = fit[1]
                beta = fit[2:(length(fit)-1)]
                sigma = fit[length(fit)]
                sigma2 = sigma*sigma

                xbeta = x_e%*%beta
                zl = (l-alpha-xbeta)/sigma
                zu = (u-alpha-xbeta)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = nr + sum(a)
                I11_22 = t(x_r)%*%x_r + t(x_e)%*%(diag(a[,1])%*%x_e)

                I11_21 = t(t(colSums(x_r))) + t(t(a)%*%x_e)
                I11_12 = t(I11_21)

                I11_33 = 2*nr + sum(c)
                I11_31 = sum(b)
                I11_13 = sum(b)
                I11_23 = t(t(b)%*%x_e)
                I11_32 = t(b)%*%x_e

                I11 = cbind(rbind(I11_11,I11_21,I11_31),
                            rbind(I11_12,I11_22,I11_32),
                            rbind(I11_13,I11_23,I11_33))

                I22 = t(g_r)%*%g_r +  t(g_e)%*%(diag(a[,1])%*%g_e)

                I21_1 = t(t(colSums(g_r))) + t(t(a)%*%g_e)
                I21_2 = t(t(x_r)%*%g_r) +  t(t(x_e)%*%(diag(a[,1])%*%g_e))
                I21_3 = t(t(b)%*%g_e)
                I21 = cbind(I21_1,I21_2,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = (t(y_r-alpha - x_r%*%beta)%*%g_r + t(y_e-alpha - xbeta +sigma*h0)%*%g_e)/sigma2
                t = s%*%ginv(Sigma)%*%t(s)
                pvalue[1,1] = pchisq(t,ng,lower.tail=FALSE)
                statistic[1,1] = t
                parameter[1,1] = ng
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }else{
                ###############################################################
                # No covariates present in the null model
                ###############################################################

                fit = epsonlyloglikmax(modeldata,c(l,u),randomindex) # Fit under H0
                alpha = fit[1]
                sigma = fit[length(fit)]
                sigma2 = sigma*sigma

                zl = (l-alpha)/sigma
                zu = (u-alpha)/sigma

                h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
                h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
                h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
                h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
                    (1-pnorm(zu)+pnorm(zl))

                a = 1 - h1 - h0*h0
                b = h0 - h2 - h0*h1
                c = 2 + 2*h1 - h3 - h1*h1

                I11_11 = nr + ne*a

                I11_33 = 2*nr + ne*c
                I11_31 = ne*b
                I11_13 = ne*b

                I11 = cbind(rbind(I11_11,I11_31),
                            rbind(I11_13,I11_33))

                I22 = t(g_r)%*%g_r +  a*t(g_e)%*%g_e

                I21_1 = t(t(colSums(g_r))) + a*t(t(colSums(g_e)))
                I21_3 = b*t(t(colSums(g_e)))
                I21 = cbind(I21_1,I21_3)
                I12 = t(I21)

                Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
                s = (t(y_r-alpha)%*%g_r + t(y_e-alpha +sigma*h0)%*%g_e)/sigma2
                t = s%*%ginv(Sigma)%*%t(s)
                pvalue[1,1] = pchisq(t,ng,lower.tail=FALSE)
                statistic[1,1] = t
                parameter[1,1] = ng
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
        }
    }
}
