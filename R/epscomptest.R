#' @title Score test EPS-complete
#'
#' @description
#' \code{epscomp.test} performs a score test for genetic variables in
#' the EPS-complete design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param SNP a matrix of variables to be tested against the null
#' @param onebyone \code{TRUE} if genetic variables should be tested one by one
#' for inclusion in the model, default \code{TRUE}
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @return Maximum likelihood estimates of the model parameters,
#' with 95 percent confidence intervals
#' \item{coef}{a vector of maximum likelihood estimates of the coefficients}
#' \item{ci}{a matrix of confidence intervals for the coefficients}
#' \item{sigma}{the estimated standard deviation}
#' \item{maf}{the estimated minor allele frequencies, returned if
#' \code{HWE = TRUE} both MAFs are unknown}
#' @details
#' The formula object is similar to that of the \code{lm} function,
#' and describes a regression model, assuming a normal distribution for
#' the residuals. See Examples.
#'
#' The data set must consist of observations of the response \code{y} of the
#' linear regression model, and (optional) any non-genetic covariates. The
#' genetic covariates (SNPs) must be observed only for the extreme-phenotype
#' individuals, and missing values for the non-extremes should be coded as
#' \code{NA}.
#'
#' If confounder = TRUE, the genetic variables are assumed to be
#' multinomally distributed, with different distribution
#' for different levels of other (non-genetic) covariates, these can
#' be specified by a vector of names \code{cx}.
#'
#' If confounder = FALSE, the distribution of \code{xg} is defined
#' by minor allele frequency (MAF)
#' and the genetic effect model.
#' @import MASS stats
#' @export
#' @examples
#' ## Create dataset:
#' N = 2000
#' xe = rnorm(n = N, mean = 2, sd = 1)
#' maf = 0.2
#' xg = sample(c(0,1,2),N,c((1-maf)^2,2*maf*(1-maf),maf^2), replace = TRUE)
#' maf2 = 0.4
#' xg2 = sample(c(0,1,2),N,c((1-maf2)^2,2*maf2*(1-maf2),maf2^2), replace = TRUE)
#' a = 50; be = 5; bg = 0.3; sigma = 2
#' y = rnorm(N, mean = a + be*xe + bg*xg, sd = sigma)
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' xg[!extreme] = NA
#' xg2[!extreme] = NA
#' ## Perform score test:
#' epscomp.test(y~xe,cbind(xg,xg2))
#' epscomp.test(y~xe,cbind(xg,xg2),onebyone = FALSE)
#'
epscomp.test = function(nullmodel, SNP, onebyone = TRUE,
                        confounder = FALSE, cx){
    if(is.null(colnames(SNP))){
        SNP = as.matrix(SNP)
        colnames(SNP) = paste0("SNP",1:dim(SNP)[2])}

    options(na.action="na.pass")
    epsdata0 = model.frame(nullmodel)
    covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    model0 = attr(terms(nullmodel),"term.labels")

    if(confounder & !onebyone){stop("The method is not available for
                                    confounding and testing multiple
                                    genetic markers simultaneously")}

    y = epsdata0[,1]
    n = length(y)

    modelnames = attr(terms(nullmodel), "term.labels")
    covariateorder = attr(terms(nullmodel), "order")
    if(length(modelnames)!= dim(covariates0)[2]){
        tonullmodel = c()
        for(i in 1:length(modelnames)){
            if(covariateorder[i] > 1){
                tonullmodel[length(tonullmodel)+1] = modelnames[i]
            }else{
                mat = as.matrix(get(all.vars(nullmodel)[(i+1)]))
                if(dim(mat)[2]>1){
                    for(j in 1:dim(mat)[2]){
                        assign(colnames(mat)[j],mat[,j])
                        tonullmodel[length(tonullmodel)+1] = colnames(mat)[j]
                    }
                }else{
                    tonullmodel[length(tonullmodel)+1] = modelnames[i]
                }
            }

        }
        # then there is a covariate in the fomula that is a matrix
        tonullmodel = unique(tonullmodel)
        nullmodel = as.formula(paste("y ~ ", paste(tonullmodel, collapse= "+")))
        options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = model.matrix(nullmodel)[,-1]
        options(na.action="na.omit")
        modelnames = attr(terms(nullmodel), "term.labels")
    }

    isxe = FALSE
    if(dim(covariates0)[2]>0){
        isxe = TRUE
        xe = covariates0
        # Confounders
        if(confounder){
            if(missing(cx)){cx = model0}
            print(paste("Confounders: ", toString(cx),sep=""))
            cind = match(cx,model0)
            xecind = as.matrix(covariates0[,cind])
            for(j in 1:length(cx)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders. \n
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

                #                 print(unique(eg_ic))
                #
                #                 for(k in 1:length(y_ic)){
                #                     if(!is.na(mean(g[x_cc[,cind] == x_ic[k,cind]])) &
                #                        !is.na(var(g[x_cc[,cind] == x_ic[k,cind]])) ){
                #                         eg_ic[k] = mean(g[x_cc[,cind] == x_ic[k,cind]])
                #                         varg_ic[k] = var(g[x_cc[,cind] == x_ic[k,cind]])
                #                         egg_ic[k] = varg_ic[k] + eg_ic[k]^2
                #                     }else{
                #                         eg_ic[k] = mean(g)
                #                         varg_ic[k] = var(g)
                #                         egg_ic[k] = varg_ic[k] + eg_ic[k]^2
                #                     }
                #                 }
                #
                #                 print(unique(eg_ic))

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
                pval = pchisq(t,1,lower.tail=F)

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
                pval = pchisq(t,1,lower.tail=F)

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
            pval = pchisq(t,ng,lower.tail=F)

            pval = pchisq(t,ng,lower.tail=F)
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
            pval = pchisq(t,ng,lower.tail=F)

            pval = pchisq(t,ng,lower.tail=F)
            result = list(t,ng,pval)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }
    }
}
