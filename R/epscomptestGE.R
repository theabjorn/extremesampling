#' @title Test for gene-environment associations under the EPS-complete design
#' @description
#' \code{epscomp.testGE} performs a likelihood ratio test for gene-environment
#' interaction variables under the EPS-complete design
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
#' @return \code{epscomp.testGE} returns
#' \item{statistic}{the value of the score test statistic}
#' \item{parameter}{the degrees of freedom of the statistic}
#' \item{p.value}{the P-value for the test}
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe+xg, which describes a regression model, y=a+be*xe+bg*xg+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The covariate xg is a SNP (single-nucleotide polymorphism).
#' The variables are taken from the environment that the
#' function is called from.
#' Both xe and xg can be matrices.
#'
#' The test considers a regression model y=a+be*xe+bg*xg+b*xe*xg+e,
#' where b=0 under the null hypothesis. The output
#' of the function gives the test statistic and p-value for the test of
#' H0: b=0. The specific gene-environment interactions that should be
#' tested is specified in \code{GE}.
#'
#' The EPS-complete design is such that the SNP genotype is only observed
#' for individuals with high and low values of the phenotype \code{y}.
#' For remaining individuals, the unobserved genotype most be coded as NA.
#' A SNP is assumed to have possible genotype 0, 1 or 2 according to the
#' number of minor-alleles. The distribution of the genotype is assumed
#' unknown and multinomially distributed. I.e. P(xg=0) = p0, P(xg=1) = p1,
#' and P(xg=2) = p2 = 1-p0-p1.
#' Hardy-Weinberg equilibrium with known or uknown MAF can be assumed,
#' then p0 = (1-q)^2, p1 = 2q(1-q) and p2 = q^2, where q is the MAF.
#' The parameters p0, p1, p2 can be given in \code{gfreq} if they are known.
#'
#' If confounder = TRUE, the genetic variables are assumed to be
#' multinomally distributed, with different distribution
#' for different levels of other (non-genetic) covariates, these can
#' be specified by a vector of names \code{cx}.
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
#' # Create the EPS-complete data set by setting
#' # the SNP values of non-extremes to NA
#' xg1[!extreme] = NA
#' xg2[!extreme] = NA
#' xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#' epscomp.testGE(y~xe1+xe2+xg1+xg2,GE = c("xe1:xg1"))$p.value

epscomp.testGE = function(nullmodel, GE, onebyone = TRUE,
                          confounder = FALSE, cx, hwe = FALSE, maf,
                          gfreq){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    covnames = colnames(covariates0)

    # which main effects are SNPs and which are environment?
    snpid = c()
    xid = c()
    ncov = length(covnames)
    for(c in 1:ncov){
        if(sum(is.na(covariates0[,c])) > 0){
            snpid[length(snpid)+1] = c
        }else if(sum(is.na(covariates0[,c])) == 0){
            xid[length(xid)+1] = c
        }
    }

    # which main effects should interact?
    interactind = list()
    nge = length(GE)
    xge = matrix(NA,ncol = nge,nrow = n)
    for(i in 1:nge){
        t = match(strsplit(GE[i],":")[[1]], covnames)
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

    if(confounder){stop("Confounding currently not allowed for interactions")}

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

