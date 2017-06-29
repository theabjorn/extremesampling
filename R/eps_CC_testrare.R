#' @title Score test EPS-CC rare variants
#' @description
#' \code{epsCC.rv.test} performs a score test for a burden of
#' rare genetic variants under the EPS complete-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param RV a matrix of genetic variants to be tested against the null
#' @param cutoffs a vector \code{c(l,u)} of the lower and upper cut-offs used
#' for extreme sampling
#' @param randomindex a binary vector that indicates if samples are random
#' or extreme
#' @param method testing the burden using \code{naive}, \code{collapsing} or
#' \code{lmm} method, see details
#' @param weights optional weights for the \code{collapsing} or
#' \code{lmm} method
#'
#' @return \code{epsCC.rv.test} returns for the whole burden of variants:
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
#' The covariate \code{xg} is a burden of several rare genetic variants.
#'
#' For the EPS complete-case design, the data is only available
#' for individuals with high and low values of the phenotype \code{y};
#' (\code{y < l} or \code{y > u}), and potentialy some randomly sampled
#' individuals. The cut-offs \code{l} and \code{u} that specify the
#' sampling must be given in the \code{cutoffs} argument.
#'
#' The \code{naive} method uses a standard score test to test the burden,
#' the \code{collaps} method tests the (weighted) sum of all variants in the
#' burden, while the \code{lmm} method is a (weighted) variance
#' component score test.
#'
#' @import MASS stats
#' @export
#' @examples
#' N = 5000 # Number of individuals in a population
#' xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
#' xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
#' maf1 = 0.01; maf2 = 0.02; maf3 = 0.005
#' xg1 = sample(c(0,1,2),N,c((1-maf1)^2,2*maf1*(1-maf1),maf1^2), replace = TRUE) # RV
#' xg2 = sample(c(0,1,2),N,c((1-maf2)^2,2*maf2*(1-maf2),maf2^2), replace = TRUE) # RV
#' xg3 = sample(c(0,1,2),N,c((1-maf3)^2,2*maf3*(1-maf3),maf3^2), replace = TRUE) # RV
#' # Model parameters
#' a = 50; be1 = 5; be2 = 8; bg1 = -0.3; bg2 = 0.6; sigma = 2
#' # Generate response y
#' y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 + bg2*xg2, sd = sigma)
#' # Identify extremes, here upper and lower 25% of population
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' # Create the EPS complete-case data set
#' y = y[extreme]; xe1 = xe1[extreme]; xe2 = xe2[extreme]
#' xg1 = xg1[extreme]; xg2 = xg2[extreme]; xg3 = xg3[extreme]
#' xg = as.matrix(cbind(xg1,xg2,xg3)); xe = as.matrix(cbind(xe1,xe2))
#'
#' # Testing
#' # Naive test
#' epsCC.rv.test(y~xe,RV=xg,cutoffs = c(l,u))
#' # Collapsing test
#' epsCC.rv.test(y~xe,RV=xg,cutoffs = c(l,u),method = "collaps",weights)
#' # Variance component test (assume linear mixed model)
#' epsCC.rv.test(y~xe,RV=xg,cutoffs = c(l,u),method = "lmm",weights)
#'

epsCC.rv.test = function(nullmodel,RV,cutoffs,randomindex,method = "naive",weights){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(length(cutoffs) != 2){stop("Invalid cutoffs vector given")}

    rsample = TRUE
    if(missing(randomindex)){
        rsample = FALSE
    }else if(sum(randomindex)==0){
        rsample = FALSE
    }

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    isx = (dim(covariates0)[2]>0) # there are covariates in the null model

    if(sum(is.na(epsdata0) > 1)){
        stop("NA values in the data not allowed")}

    l = min(cutoffs)
    u = max(cutoffs)

    RV = as.matrix(RV)
    # if(dim(RV)[2]<2){
    #     stop("RV must be a set of more than one genetic variant")
    # }

    if(missing(weights)){
        weights = rep(1,dim(RV)[2])
    }

    if(method == "naive"){
        epsCC.rv.test.naive(epsdata0,covariates0,RV,isx,l,u,rsample, randomindex)
    }else if(method == "collaps"){
        g = rep(0,dim(RV)[1])
        for(c in 1:dim(RV)[2]){
            g = g + c(weights[c])*c(RV[,c])
        }
        epsCC.test(nullmodel,SNP=g,cutoffs,randomindex)
    }else if(method == "lmm"){
        epsCC.rv.test.lmm(epsdata0,covariates0,RV,isx,l,u,rsample,randomindex,weights)
    }
}

