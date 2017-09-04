#' @title Score test EPS-AC rare variants
#' @description
#' \code{epsAC.rv.test} performs a score test for a burden of
#' rare genetic variants under the EPS all-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param xg a matrix of variables to be tested against the null (NA for
#' not genotyped individuals)
#' @param confounder (optional) vector of names of confounding (non-genetic) covariates
#' @param method testing the burden using \code{naive}, \code{collapsing} or
#' \code{varcomp} method, see details
#' @param weights optional weights for the \code{collapsing} or
#' \code{varcomp} method
#'
#' @return \code{epsAC.rv.test} returns for the whole burden of variants:
#' \item{statistic}{the score test statistic}
#' \item{p.value}{the P-value}
#'
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe, which describes a regression model, y=a+be*xe+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The variables are taken from the environment that the
#' function is called from.
#' The null hypothesis bg=0 is tested for the model y=a+be*xe+bg*xg+e.
#' The covariate \code{xg} is a burden of several rare genetic variants,
#' where missing values are coded as NA. Missing-mechanism must be MCAR or MAR.
#' Missing-mechanism assumed equal for all genetic variants.
#'
#' Confounders are discrete covariates (xe) and the distribution of xg is
#' modelled for each level of unique value of xe.
#'
#' The \code{naive} method uses a standard score test to test the burden,
#' the \code{collaps} method tests the (weighted) sum of all variants in the
#' burden, while the \code{varcomp} method is a (weighted) variance
#' component score test.
#'
#' @import MASS stats
#' @export
#' @examples
#' N = 5000 # Number of individuals in a population
#' xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
#' xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
#' maf1 = 0.01; maf2 = 0.02; maf3 = 0.005
#' xg1 = sample(c(0,1,2),N,c((1-maf1)^2,2*maf1*(1-maf1),maf1^2), replace = TRUE) # xg
#' xg2 = sample(c(0,1,2),N,c((1-maf2)^2,2*maf2*(1-maf2),maf2^2), replace = TRUE) # xg
#' xg3 = sample(c(0,1,2),N,c((1-maf3)^2,2*maf3*(1-maf3),maf3^2), replace = TRUE) # xg
#' # Model parameters
#' a = 50; be1 = 5; be2 = 8; bg1 = -0.3; bg2 = 0.6; sigma = 2
#' # Generate response y
#' y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 + bg2*xg2, sd = sigma)
#' # Identify extremes, here upper and lower 25% of population
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' # Create the EPS all-case data set
#' xg1[!extreme] = NA; xg2[!extreme] = NA; xg3[!extreme] = NA;
#' xg = as.matrix(cbind(xg1,xg2,xg3))
#' xe = as.matrix(cbind(xe1,xe2))
#'
#' # Testing
#' # Naive test
#' epsAC.rv.test(y~xe,xg=xg)
#' # Collapsing test
#' epsAC.rv.test(y~xe,xg=xg,method = "collaps",weights)
#' # Variance component test (assume linear mixed model)
#' epsAC.rv.test(y~xe,xg=xg,method = "varcomp",weights)
#'

epsAC.rv.test = function(nullmodel,xg,confounder,method = "naive",weights){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    conf = FALSE
    if(!missing(confounder)){conf = TRUE}

    isxe = FALSE
    if(dim(covariates0)[2]>0){
        isxe = TRUE
        xe = covariates0
        if(dim(covariates0)[2] ==1){
            colnames(covariates0) = colnames(epsdata0)[2]
        }
        # Confounders
        if(conf){
            if(is.na(match(confounder,colnames(covariates0)))){
                stop("The name of the confounder given is not the name of a covariate.")
            }
            message(paste("Confounders: ", toString(confounder),sep=""))
            cind = match(confounder,colnames(covariates0))
            xecind = as.matrix(covariates0[,cind])
            for(j in 1:length(confounder)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                         Please recode you confounder to satisfy this.")
                }
                }
            }
        }

    xg = as.matrix(xg)
    #if(dim(xg)[2]<2){
    #    stop("xg must be a set of more than one genetic variant")
    #}

    if(missing(weights)){
        weights = rep(1,dim(xg)[2])
    }

    if(method == "naive"){
        epsAC.rv.test.naive(epsdata0,covariates0,xg,isxe,conf,cind)
    }else if(method == "collaps"){
        g = rep(0,dim(xg)[1])
        for(c in 1:dim(xg)[2]){
            g = g + c(weights[c])*c(xg[,c])
        }
        epsAC.test(nullmodel,xg=g,confounder)
    }else if(method == "varcomp"){
        epsAC.rv.test.varcomp(epsdata0,covariates0,xg,isxe,conf,cind,weights)
    }
}

