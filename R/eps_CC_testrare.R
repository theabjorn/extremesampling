#' @title Score test EPS-CC rare variants
#' @description
#' \code{epsCC.rv.test} performs a score test for a burden of
#' rare genetic variants under the EPS complete-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param xg a matrix of genetic variants to be tested against the null
#' @param l cutoff for lower extreme, can be sample-specific or general
#' @param u cutoff for upper extreme, can be sample-specific or general
#' @param method testing the burden using \code{simple}, \code{collapsing} or
#' \code{varcomp} method, see details
#' @param weights optional weights
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
#' The \code{simple} method uses a standard score test to test the burden,
#' the \code{collapse} method tests the (weighted) sum of all variants in the
#' burden, while the \code{varcomp} method is a (weighted) variance
#' component score test.
#'
#' The \code{varcomp} method is a special case of the popular SKAT method
#' for a linear weighted kernel under extreme phenotype sampling. The p-value is found
#' using the function \code{\link[CompQuadForm]{davies}} in the CompQuadForm
#' package.
#'
#' @references
#' \insertRef{quadRpackage}{extremesampling},
#'
#' \insertRef{wu2011SKAT}{extremesampling}
#'
#' @seealso \code{\link[SKAT]{SKAT}} for the SKAT test for random samples,
#' \code{\link[CompQuadForm]{davies}} for the Davies method,
#' \code{\link{epsCC.test}} for a common variant SNP by SNP test for
#' complete-case extreme sampling
#'
#' @import MASS stats CompQuadForm
#' @export
#' @examples
#'
#' N = 1000
#' # Generate environmental covariates
#' xe1 = rbinom(N,1,0.5); xe2 = rnorm(N,2,1)
#' # Generate genetic covariates (common variants)
#' cv1 = rbinom(N,2,0.2); cv2 = rbinom(N,2,0.2)
#' # Generate phenotype
#' y = rnorm(N, mean = 1 + 2*xe1 + 3*xe2 + 0.5*cv1 + 0.1*cv2,2)
#' # Define extremes
#' u = quantile(y,probs = 3/4,na.rm=TRUE); l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' cv1[!extreme] = NA; cv2[!extreme] = NA;
#'
#' # Complete case data set
#' xe_CC = cbind(xe1[extreme], xe2[extreme])
#' xg_CC = cbind(cv1[extreme], cv2[extreme])
#' y_CC = y[extreme]
#'
# # Simple test
#' epsCC.rv.test(y_CC ~ xe_CC,xg = xg_CC,l,u)
#' # Collapsing test
#' epsCC.rv.test(y_CC ~ xe_CC,xg = xg_CC,method = "collapse",l,u)
#' # Variance component test
#' epsCC.rv.test(y_CC ~ xe_CC,xg = xg_CC,method = "varcomp",l,u)
#'

epsCC.rv.test = function(nullmodel,xg,l,u,method = "simple",weights){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        xe = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    isx = (dim(xe)[2]>0) # there are covariates in the null model

    if(sum(is.na(epsdata0) > 1)){
        stop("NA values in the data not allowed")}

    y = epsdata0[,1]

    xg = as.matrix(xg)
    # if(dim(xg)[2]<2){
    #     stop("xg must be a set of more than one genetic variant")
    # }

    if(missing(weights)){
        message("weights = 1")
        weights = rep(1,dim(xg)[2])
    }

    if(isx){
        if(method == "simple"){
            eps_CC_test_simple_EX(y,xe,xg,l,u)
        }else if(method == "collapse"){
            g = colSums(t(xg)*(weights)^2)
            epsCC.test(y~xe,xg=g,l,u)
        }else if(method == "varcomp"){
            wxg = t(t(xg)*weights)
            eps_CC_test_varcomp_EX(y,xe,wxg,l,u)
        }
    }else{
        if(method == "simple"){
            eps_CC_test_simple_E(y,xg,l,u)
        }else if(method == "collapse"){
            g = colSums(t(xg)*(weights)^2)
            epsCC.test(nullmodel,xg=g,l,u)
        }else if(method == "varcomp"){
            wxg = t(t(xg)*weights)
            eps_CC_test_varcomp_E(y,wxg,l,u)
        }
    }



}

