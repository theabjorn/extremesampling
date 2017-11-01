#' @title Score test EPS-AC rare variants
#' @description
#' \code{epsAC.rv.test} performs a score test for a burden of
#' rare genetic variants under the EPS all-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param xg a matrix of variables to be tested against the null (NA for
#' not genotyped individuals)
#' @param confounder (optional) vector of names of confounding (non-genetic) covariates
#' @param method testing the burden using \code{simple}, \code{collapse} or
#' \code{varcomp} method, see details
#' @param weights optional weights
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
#' \code{\link{epsAC.test}} for a common variant SNP by SNP test for
#' all-case extreme sampling
#'
#' @import MASS stats CompQuadForm
#' @export
#' @examples
#' # Simple test
#' epsAC.rv.test(phenoRV~EPSxe,xg=gRV)
#' # Collapsing test
#' epsAC.rv.test(phenoRV~EPSxe,xg=gRV,method = "collapse")
#' # Variance component test
#' epsAC.rv.test(phenoRV~EPSxe,xg=gRV,method = "varcomp")
#'
#' Control for include confounders
#' epsAC.rv.test(phenoRV~EPSxe1 + EPSxe2,xg=gRV,confounder = EPSxe1)

epsAC.rv.test = function(nullmodel,xg,confounder,method = "simple",weights){
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
    ng = dim(xg)[2]

    conf = FALSE
    if(!missing(confounder)){
        conf = TRUE
        xec = as.matrix(confounder)
        for(j in 1:dim(xec)[2]){
            if(length(unique(xec[,j])) > 10){
                stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
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

    if(missing(weights)){
        weights = rep(1,dim(xg)[2])
    }

    if(isxe){
        if(!conf){
            if(method == "simple"){
                res = eps_AC_test_rv_X(y,xe,xg,weights)
                t = crossprod(res[[1]],crossprod(ginv(res[[2]]),res[[1]]))
                pvalue = pchisq(t,ng,lower.tail=FALSE)
                statistic = t
                result = list(statistic,pvalue)
                names(result) = c("statistic","p.value")
                return(result)
            }else if(method == "collapse"){
                g = colSums(t(xg)*(weights)^2)
                epsAC.test(y~xe,xg=g)
            }else if(method == "varcomp"){
                res = eps_AC_test_rv_X(y,xe,xg,weights)
                d = eigen(res[[2]])$values
                teststat = crossprod(res[[1]])
                pvalue = davies(teststat,d,acc = 0.00005)$Qq
                result = list(teststat,pvalue,d)
                names(result) = c("statistic","p.value","eigen")
                return(result)
            }
        }else{
            if(method == "simple"){
                message(paste("Confounding assumed"))
                res = eps_AC_test_rv_X_conf(y,xe,xec,xg,weights)
                t = crossprod(res[[1]],crossprod(ginv(res[[2]]),res[[1]]))
                pvalue = pchisq(t,ng,lower.tail=FALSE)
                statistic = t
                result = list(statistic,pvalue)
                names(result) = c("statistic","p.value")
                return(result)
            }else if(method == "collapse"){
                g = colSums(t(xg)*(weights)^2)
                epsAC.test(y~xe,xg=g,confounder = confounder)
            }else if(method == "varcomp"){
                message(paste("Confounding assumed"))
                res = eps_AC_test_rv_X_conf(y,xe,xec,xg,weights)
                d = eigen(res[[2]])$values
                teststat = crossprod(res[[1]])
                pvalue = davies(teststat,d,acc = 0.00005)$Qq
                result = list(teststat,pvalue,d)
                names(result) = c("statistic","p.value","eigen")
                return(result)
            }
        }
    }else{
        if(method == "simple"){
            res = eps_AC_test_rv_noX(y,xg,weights)
            t = crossprod(res[[1]],crossprod(ginv(res[[2]]),res[[1]]))
            pvalue = pchisq(t,ng,lower.tail=FALSE)
            statistic = t
            result = list(statistic,pvalue)
            names(result) = c("statistic","p.value")
            return(result)
        }else if(method == "collapse"){
            g =colSums(t(xg)*(weights)^2)
            epsAC.test(y~1,xg=g)
        }else if(method == "varcomp"){
            res = eps_AC_test_rv_noX(y,xg,weights)
            d = eigen(res[[2]])$values
            teststat = crossprod(res[[1]])
            pvalue = davies(teststat,d,acc = 0.00005)$Qq
            result = list(teststat,pvalue,d)
            names(result) = c("statistic","p.value","eigen")
            return(result)
        }
    }
}

