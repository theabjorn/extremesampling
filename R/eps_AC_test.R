#' @title Score test EPS-AC
#' @description
#' \code{epsAC.test} performs a score test for common genetic variants
#' under the EPS all-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis, see details
#' @param xg a matrix of variables to be tested against the null (NA for
#' not genotyped individuals)
#' @param confounder confounding (non-genetic) covariates
#' @return \code{epsAC.test} returns
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
#' The covariate xg is one or more, genetic markers where missing values
#' are coded as NA. Missing-mechanism must be MCAR or MAR.
#' If xg is a matrix, each variant (column) is tested against the null model.
#'
#' Confounders are discrete covariates (xe) and the distribution of xg is
#' modelled for each level of unique value of xe.
#'
#' @import MASS stats
#' @export
#' @examples
#' # Testing
#' epsAC.test(phenoCV ~ EPSxe,xg = gCV)
#' epsAC.test(phenoCV ~ EPSxe1 + EPSxe2,xg = gCV,confounder = EPSxe1)

epsAC.test = function(nullmodel, xg, confounder){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(is.null(colnames(xg))){
        xg = as.matrix(xg)
        colnames(xg) = paste0("xg",1:dim(xg)[2])}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        xe = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    conf = FALSE

    if(!missing(confounder)){
        conf = TRUE
        xec = as.matrix(confounder)
        if(dim(xec)[1]!=n){
            stop("Confounder has incorrect dimension")
        }
        for(j in 1:dim(xec)[2]){
        if(length(unique(xec[,j])) > 50){
            stop("Only discrete confounders with less than or equal to 50 unique levels are accepted as confounders.
                 Please recode you confounder to satisfy this.")
            }
        }
    }

    isxe = FALSE
    if(dim(xe)[2]>0){
        isxe = TRUE
        if(is.null(colnames(xe))){
            colnames(xe) = paste0("xe",1:dim(xe)[2])}
    }


    if(isxe & !conf){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # No confounding assumed
        ###############################################################

        eps_AC_test_x(y,xe,xg)

    }else if(isxe & conf){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # Confounding assumed
        ###############################################################
        eps_AC_test_x_conf(y,xe,xec,xg)

    }else{
        ##############################################################
        # no environmental covariates (xe) present in the null model
        ##############################################################

        eps_AC_test_nox(y,xg)

    }

}
