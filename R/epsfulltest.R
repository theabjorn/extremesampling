#' @title Test for associations under the EPS-full design
#' @description
#' \code{epsfull.test} performs a score test for genetic variables in
#' the EPS-full design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis, see details
#' @param SNP a matrix of variables to be tested against the null
#' @param onebyone \code{TRUE} if genetic variables should be tested one by one
#' for inclusion in the model, default \code{TRUE}
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @return \code{epsfull.test} returns
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
#' The EPS-full design is such that the SNP genotype is only observed
#' for individuals with high and low values of the phenotype \code{y}.
#' For remaining individuals, the unobserved genotype most be coded as NA.
#' A SNP is assumed to have possible genotype 0, 1 or 2 according to the
#' number of minor-alleles.
#'
#' If confounder = TRUE, the genetic variables are assumed to have
#' different distribution for different levels of other (non-genetic)
#' covariates, these can be specified by a vector of names \code{cx}.
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
#' # Create the EPS-full data set by setting
#' # the SNP values of non-extremes to NA
#' xg1[!extreme] = NA
#' xg2[!extreme] = NA
#' xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#' epsfull.test(y~xe1+xe2,SNP = xg)$p.value
#' epsfull.test(y~xe,SNP = xg,onebyone=FALSE)$p.value

epsfull.test = function(nullmodel, SNP, onebyone = TRUE,
                        confounder = FALSE, cx){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(is.null(colnames(SNP))){
        SNP = as.matrix(SNP)
        colnames(SNP) = paste0("SNP",1:dim(SNP)[2])}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    if(confounder & !onebyone){stop("The method is not available for
                                    confounding and testing multiple
                                     genetic markers simultaneously")}

    y = epsdata0[,1]
    n = length(y)

    isxe = FALSE
    if(dim(covariates0)[2]>0){
        isxe = TRUE
        xe = covariates0
        # Confounders
        if(confounder){
            if(missing(cx)){cx = colnames(covariates0)
            }else if(is.na(match(cx,colnames(covariates0)))){
                stop("The name of the confounder given is not the name of a covariate.")
            }
            message(paste("Confounders: ", toString(cx),sep=""))
            cind = match(cx,colnames(covariates0))
            xecind = as.matrix(covariates0[,cind])
            for(j in 1:length(cx)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                         Please recode you confounder to satisfy this.")
                }
            }
        }
    }

    totest = colnames(SNP)
    xg = as.matrix(SNP)
    ng = dim(xg)[2]
    extreme = !is.na(xg[,1])
    y_cc = y[extreme]
    y_ic = y[!extreme]
    g_cc = as.matrix(xg[extreme,])
    if(isxe){
        x_cc = as.matrix(xe[extreme,])
        x_ic = as.matrix(xe[!extreme,])
    }
    n_cc = length(y_cc)
    n_ic = length(y_ic)

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

        if(isxe & !confounder){
            ###############################################################
            # Environmental covariates (xe) present in the null model
            # No confounding assumed
            ###############################################################
            fit = lm(y~xe) # Fit under H0
            coefs = fit$coef
            alpha = coefs[1]
            beta = coefs[2:length(coefs)]
            sigma = summary(fit)$sigma
            sigma2 = sigma*sigma

            for(i in 1:ng){
                extreme = !is.na(xg[,i])
                y_cc = y[extreme]
                y_ic = y[!extreme]
                g = as.matrix(xg[extreme,i])
                x_cc = as.matrix(xe[extreme,])
                x_ic = as.matrix(xe[!extreme,])
                n_cc = length(y_cc)
                n_ic = length(y_ic)

                eg = mean(g)
                varg = var(g)
                egg = mean(g*g)

                I11 = cbind(c(n,colSums(xe),2*sum(y-alpha-xe%*%beta)/sigma),
                            rbind(colSums(xe),t(xe)%*%xe,
                                  2*(t(y-alpha-xe%*%beta)%*%xe)/sigma),
                            rbind(2*sum(y-alpha-xe%*%beta)/sigma,
                                  2*t(xe)%*%(y-alpha-xe%*%beta)/sigma,
                                  ((3/sigma2)*sum((y-alpha-xe%*%beta)^2)-n)))

                I22 = sum(g^2) + n_ic*egg -
                    (1/sigma2)*varg*sum((y_ic-alpha-x_ic%*%beta)^2)

                I12_1 = n_cc*mean(g) + n_ic*eg

                I12_2 = t(x_cc)%*%g + colSums(x_ic)*eg

                I12_3 = 2*sum((y_cc-alpha-x_cc%*%beta)*g/sigma) +
                    2*sum((y_ic-alpha-x_ic%*%beta)*eg/sigma)

                I12 = matrix(c(I12_1,I12_2,I12_3))

                Sigma = (1/sigma2)*(I22 - t(I12)%*%ginv(I11)%*%I12)
                s = (sum((y_cc - alpha - x_cc%*%beta)*g) +
                         sum((y_ic - alpha-x_ic%*%beta)*eg))/sigma2
                t = (s*s)/Sigma
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[i,] = t
                parameter[i,] = 1
                pvalue[i,] = pval
            }
            result = list(statistic,parameter,pvalue)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }else if(isxe & confounder){
            ###############################################################
            # Environmental covariates (xe) present in the null model
            # Confounding assumed
            ###############################################################
            fit = lm(y~xe) # Fit under H0
            coefs = fit$coef
            alpha = coefs[1]
            beta = coefs[2:length(coefs)]
            sigma = summary(fit)$sigma
            sigma2 = sigma*sigma

            # confounding
            ux = as.matrix(unique(x_cc[,cind]))
            nu = dim(ux)[1]
            if(nu != dim(as.matrix(unique(xe[,cind])))[1]){
                warning("All unique levels of confounder not found in extreme sample")
            }
            uindex = list()
            for(u in 1:nu){
                uindex[[u]] = (x_cc[,cind] == ux[u,])
            }

            for(i in 1:ng){
                extreme = !is.na(xg[,i])
                y_cc = y[extreme]
                y_ic = y[!extreme]
                g = as.matrix(xg[extreme,i])
                x_cc = as.matrix(xe[extreme,])
                x_ic = as.matrix(xe[!extreme,])
                n_cc = length(y_cc)
                n_ic = length(y_ic)

                eg_ic = c()
                varg_ic = c()
                egg_ic = c()

                for(u in 1:nu){
                    uind = uindex[[u]]
                    tempind = (x_ic[,cind] == ux[u,])
                    if(!is.na(mean(g[uind])) &
                       !is.na(var(g[uind]))){
                        eg_ic[tempind] = mean(g[uind])
                        varg_ic[tempind] = var(g[uind])
                        egg_ic[tempind] = varg_ic[tempind] + eg_ic[tempind]^2
                    }else{
                        eg_ic[tempind] = mean(g)
                        varg_ic[tempind] = var(g)
                        egg_ic[tempind] = varg_ic[tempind] + eg_ic[tempind]^2
                    }
                }

                I11 = cbind(c(n,colSums(xe),2*sum(y-alpha-xe%*%beta)/sigma),
                            rbind(colSums(xe),t(xe)%*%xe,
                                  2*(t(y-alpha-xe%*%beta)%*%xe)/sigma),
                            rbind(2*sum(y-alpha-xe%*%beta)/sigma,
                                  2*t(xe)%*%(y-alpha-xe%*%beta)/sigma,
                                  ((3/sigma2)*sum((y-alpha-xe%*%beta)^2)-n))
                )

                I22 = sum(g^2) + sum(egg_ic) -
                    (1/sigma2)*sum(varg_ic*(y_ic-alpha-x_ic%*%beta)^2)

                I12_1 = sum(g) + sum(eg_ic)

                I12_2 = t(x_cc)%*%g + t(x_ic)%*%eg_ic

                I12_3 = 2*sum((y_cc-alpha-x_cc%*%beta)*g/sigma) +
                    2*sum((y_ic-alpha-x_ic%*%beta)*eg_ic/sigma)

                I12 = matrix(c(I12_1,I12_2,I12_3))

                Sigma = (1/sigma2)*(I22 - t(I12)%*%ginv(I11)%*%I12)
                s = (sum((y_cc - alpha - x_cc%*%beta)*g) +
                         sum((y_ic - alpha-x_ic%*%beta)*eg_ic))/sigma2
                t = (s*s)/Sigma
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
            # no environmental covariates (xe) present in the null model
            ##############################################################
            fit = lm(y~1) # Fit under H0
            coefs = fit$coef
            alpha = mean(y)
            sigma = sd(y)
            sigma2 = sigma*sigma

            for(i in 1:ng){
                extreme = !is.na(xg[,i])
                y_cc = y[extreme]
                y_ic = y[!extreme]
                g = as.matrix(xg[extreme,i])
                n_cc = length(y_cc)
                n_ic = length(y_ic)

                eg = mean(g)
                varg = var(g)
                egg = mean(g*g)

                I11 = cbind(c(n,2*sum(y-alpha)/sigma),
                            c(2*sum(y-alpha)/sigma,
                              ((3/sigma2)*sum((y-alpha)^2)-n)))

                I22 = n_cc*mean(g*g) + n_ic*egg -
                    (n_ic/sigma2)*varg*mean((y_ic-alpha)^2)
                I12_1 = n_cc*mean(g) + n_ic*eg
                I12_3 = 2*sum((y_cc-alpha)*g/sigma) +
                    2*sum((y_ic-alpha)*eg/sigma)
                I12 = matrix(c(I12_1,I12_3))

                Sigma = (1/sigma2)*(I22 - t(I12)%*%ginv(I11)%*%I12)
                s = (sum((y_cc - alpha)*g) + sum((y_ic - alpha)*eg))/sigma2
                t = (s*s)/Sigma
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
        ###############################################################
        # Test genetic variables all together
        ###############################################################

        if(isxe){
            ##############################################################
            # Environmental covariates (xe) present in the null model
            ##############################################################
            fit = lm(y~xe) # fit under H0
            coefs = fit$coef
            alpha = coefs[1]
            beta = coefs[2:length(coefs)]
            sigma = summary(fit)$sigma
            sigma2 = sigma*sigma

            eg = matrix(colMeans(g_cc),nrow = 1)
            egg = (1/n_cc)*t(g_cc)%*%g_cc
            varg = var(g_cc)

            I11 = cbind(c(n,colSums(xe),2*sum(y-alpha-xe%*%beta)/sigma),
                        rbind(colSums(xe),t(xe)%*%xe,
                              2*(t(y-alpha-xe%*%beta)%*%xe)/sigma),
                        rbind(2*sum(y-alpha-xe%*%beta)/sigma,
                              2*t(xe)%*%(y-alpha-xe%*%beta)/sigma,
                              ((3/sigma2)*sum((y-alpha-xe%*%beta)^2)-n)))

            I22 = t(g_cc)%*%g_cc + n_ic*egg -
                (n_ic/sigma2)*varg*mean((y_ic-alpha-x_ic%*%beta)^2)

            I12_1 = colSums(g_cc) + n_ic*eg

            I12_2 = t(x_cc)%*%g_cc + colSums(x_ic)%*%eg

            I12_3 = 2*(t(y_cc-alpha-x_cc%*%beta)%*%g_cc)/sigma +
                2*sum(y_ic-alpha-x_ic%*%beta)*eg/sigma

            I12 = rbind(I12_1,I12_2,I12_3)

            Sigma = (1/sigma2)*(I22 - t(I12)%*%ginv(I11)%*%I12)
            s = (t(y_cc - alpha - x_cc%*%beta)%*%g_cc +
                     sum(y_ic - alpha-x_ic%*%beta)*eg)/sigma2
            t = s%*%ginv(Sigma)%*%t(s)
            pval = pchisq(t,ng,lower.tail=FALSE)

            pval = pchisq(t,ng,lower.tail=FALSE)
            result = list(t,ng,pval)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }else{
            ##############################################################
            # no environmental covariates (xe) present in the null model
            ##############################################################
            fit = lm(y~1) # Fit under H0
            coefs = fit$coef
            alpha = coefs[1]
            sigma = summary(fit)$sigma
            sigma2 = sigma*sigma

            eg = matrix(colMeans(g_cc),nrow = 1)
            egg = (1/n_cc)*t(g_cc)%*%g_cc
            varg = var(g_cc)

            I11 = cbind(c(n,2*sum(y-alpha)/sigma),
                        rbind(2*sum(y-alpha)/sigma,
                              ((3/sigma2)*sum((y-alpha)^2)-n)))

            I22 = t(g_cc)%*%g_cc + n_ic*egg -
                (n_ic/sigma2)*varg*mean((y_ic-alpha)^2)

            I12_1 = colSums(g_cc) + n_ic*eg

            I12_3 = 2*(t(y_cc-alpha)%*%g_cc)/sigma + 2*sum(y_ic-alpha)*eg/sigma

            I12 = rbind(I12_1,I12_3)

            Sigma = (1/sigma2)*(I22 - t(I12)%*%ginv(I11)%*%I12)
            s = (t(y_cc - alpha)%*%g_cc + sum(y_ic - alpha)*eg)/sigma2
            t = s%*%ginv(Sigma)%*%t(s)
            pval = pchisq(t,ng,lower.tail=FALSE)

            pval = pchisq(t,ng,lower.tail=FALSE)
            result = list(t,ng,pval)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }
    }
}
