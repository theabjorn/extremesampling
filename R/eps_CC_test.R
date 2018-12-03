#' @title Score test EPS-CC
#' @description
#' \code{epsCC.test} performs a score test for common genetic variants
#' under the EPS complete-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param xg a matrix of genetic variants to be tested against the null
#' @param l cutoff for lower extreme, can be sample-specific or general
#' @param u cutoff for upper extreme, can be sample-specific or general
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
#' The covariate xg is one or more genetic markers.
#'
#' Variables are only available
#' for individuals with high and low values of the phenotype \code{y};
#' (\code{y < l} or \code{y > u}), and potentialy some randomly sampled
#' individuals. The cut-offs \code{l} and \code{u} that specify the
#' sampling must be given in the \code{cutoffs} argument.
#'
#' @import MASS stats
#' @export
#' @examples
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
#' epsCC.test(y_CC ~ xe_CC,xg = xg_CC,l,u)
#'

epsCC.test = function(nullmodel,xg,l,u){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}


    xg = as.matrix(xg)
    if(is.null(colnames(xg))){
        colnames(xg) = paste0("xg",1:dim(xg)[2])}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        xe = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]

    isx = (dim(xe)[2]>0)

    if(sum(is.na(epsdata0) > 1)){
        stop("NA values in the data not allowed")}

    ###################################################################
    # No random sample, extreme-phenotype individuals only
    ###################################################################
    if(sum((y>l & y<u))>0){
        stop("Incorrect data format")}

    if(isx){
        ###############################################################
        # Covariates present in the null model
        ###############################################################
        eps_CC_test_ex(y,xe,xg,l,u)

    }else{
        ##############################################################
        # No covariates present in the null model
        ##############################################################
        # = score test for random sample
        #stop("No covariates, use score test for random samples")
        eps_CC_test_e(y,xg,l,u)
    }
}
