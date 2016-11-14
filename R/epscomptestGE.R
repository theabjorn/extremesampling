#' @title Likelihood ratio test EPS-complete
#'
#' @description
#' \code{epscomp.lrtest} performs a likelihood ratio test for genetic
#' variables in the EPS-complete design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param GE a list of interactions, with the colon-symbol used to denote
#' interaction, e.g. \code{c("xe1:xg1","xe1:xg2")}
#' @param onebyone \code{TRUE} if genetic variables should be tested one by one
#' for inclusion in the model, default \code{TRUE}
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @param hwe \code{TRUE} if Hardy-Weinberg equilibrium is assumed, default
#' set to \code{FALSE}
#' @param maf optional value for the minor allele frequencies under HWE
#' @param gfreq frequencies of genotypes for each SNP if known or
#' otherwise estimated
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
#'
epscomp.testGE = function(nullmodel, GE, onebyone = TRUE,
                          confounder = FALSE, cx, hwe = FALSE, maf,
                          gfreq){
    options(na.action="na.pass")
    epsdata0 = model.frame(nullmodel)
    covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    modelnames = attr(terms(nullmodel), "term.labels")
    if(length(modelnames)!= dim(covariates0)[2]){
        toformula = c()
        for(i in 1:length(modelnames)){
            mat = as.matrix(get(all.vars(nullmodel)[(i+1)]))
            if(dim(mat)[2]>1){
                for(j in 1:dim(mat)[2]){
                    assign(colnames(mat)[j],mat[,j])
                    toformula[length(toformula)+1] = colnames(mat)[j]
                }
            }else{
                toformula[length(toformula)+1] = modelnames[i]
            }

        }
        # then there is a covariate in the fomula that is a matrix
        nullmodel = as.formula(paste("y ~ ", paste(toformula, collapse= "+")))
        options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = model.matrix(nullmodel)[,-1]
        options(na.action="na.omit")
        modelnames = attr(terms(nullmodel), "term.labels")
    }

    model0 = attr(terms(nullmodel),"term.labels")

    y = epsdata0[,1]
    n = length(y)
    interactind = list()

    covariateorder = attr(terms(nullmodel), "order")
    maineffects = attr(terms(nullmodel),"term.labels")[covariateorder == 1]

    snpid = c()
    xid = c()
    nmain = sum(covariateorder == 1)
    for(c in 1:nmain){
        if(sum(is.na(covariates0[,c])) > 0){
            snpid[length(snpid)+1] = c
        }else if(sum(is.na(covariates0[,c])) == 0){
            xid[length(xid)+1] = c
        }
    }

    nge = length(GE)
    xge = matrix(NA,ncol = nge,nrow = n)
    for(i in 1:nge){
        t = match(strsplit(GE[i],":")[[1]], maineffects)
        if(sum(is.na(t)>0)){
            stop(paste("Cannot include interaction term ",toString(GE[i]),
                       ", because corresponding main effects were not found in the null model",
                       sep = ""))
        }
        xge[,i] = covariates0[,t[1]]*covariates0[,t[2]]
        if(t[1] %in% snpid){
            # always genetic variant first
            interactind[[i]] = c(match(t[1],snpid),match(t[2],xid))
        }else if(t[1] %in% xid){
            interactind[[i]] = c(match(t[2],snpid),match(t[1],xid))
        }
    }

    if(confounder){stop("Confounding currently not allowed for
                        the likelihood ratio rest")}
    geneffect = "additive"
    if(missing(gfreq)){gfreq = NA}

    xg = as.matrix(covariates0[,snpid])
    xe = as.matrix(covariates0[,xid])
    data = cbind(y,xe,xg)
    ng = dim(xg)[2]

    if(onebyone){
        #######################################################################
        # Test one interaction at a time
        #######################################################################
        nint = length(interactind)

        statistic = matrix(NA,ncol = nint, nrow = 1)
        parameter = matrix(NA,ncol = nint, nrow = 1)
        pvalue = matrix(NA,ncol = nint, nrow = 1)
        colnames(statistic) = GE
        colnames(parameter) = GE
        colnames(pvalue) = GE

        if(hwe){
            #######################################################
            # Hardy Weinberg
            #######################################################
            if(missing(maf)){maf = NA}
            fit0 = epscomploglikmax(data,ng, hwe = TRUE, maf = maf,
                                    geneffect = geneffect,
                                    ll = TRUE)[[1]]
            for(l in 1:nint){
                fit1 = epscomploglikmaxint(data,ng = ng,
                                           hwe = TRUE, maf = maf,
                                           interactind = list(interactind[[l]]),
                                           ll = TRUE)[[1]]
                t = -2*(fit0 - fit1)
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[,l] = t
                parameter[,l] = 1
                pvalue[,l] = pval
            }
            result = list(statistic,parameter,pvalue)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }else{
            #######################################################
            # Hardy Weinberg not assumed
            #######################################################
            fit0 = epscomploglikmax(data,ng,
                                    geneffect = geneffect,
                                    ll = TRUE)[[1]]
            for(l in 1:nint){
                fit1 = epscomploglikmaxint(data,ng = ng,
                                           interactind = list(interactind[[l]]),
                                           ll = TRUE)[[1]]
                t = -2*(fit0 - fit1)
                pval = pchisq(t,1,lower.tail=FALSE)

                statistic[,l] = t
                parameter[,l] = 1
                pvalue[,l] = pval
            }
            result = list(statistic,parameter,pvalue)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }
    }else{
        #######################################################################
        # Test all interactions simultaneously
        #######################################################################
        if(hwe){
            #######################################################
            # Hardy Weinberg
            #######################################################
            if(missing(maf)){maf = NA}
            fit0 = epscomploglikmax(data,ng, hwe = TRUE, maf = maf,
                                    geneffect = geneffect,
                                    ll = TRUE)[[1]]

            fit1 = epscomploglikmaxint(data,ng = ng,
                                       hwe = TRUE, maf = maf,
                                       interactind = interactind,
                                       ll = TRUE)[[1]]
            t = -2*(fit0 - fit1)
            pval = pchisq(t,ng,lower.tail=FALSE)

            result = list(t,ng,pval)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }else{
            #######################################################
            # Hardy Weinberg not assumed
            #######################################################
            fit0 = epscomploglikmax(data,ng,
                                    geneffect = geneffect,
                                    ll = TRUE)[[1]]

            fit1 = epscomploglikmaxint(data,ng = ng,
                                       interactind = interactind,
                                       ll = TRUE)[[1]]
            t = -2*(fit0 - fit1)
            pval = pchisq(t,ng,lower.tail=FALSE)

            result = list(t,ng,pval)
            names(result) = c("statistic","parameter","p.value")
            return(result)
        }
    }
}
