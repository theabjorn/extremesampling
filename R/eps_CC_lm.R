#' Fit linear model to EPS-only data
#' @description
#' \code{epsCC.lm} fits a normal linear regression model to EPS-only data
#' @param formula an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model to be fitted, see details
#' @param l cutoff for lower extreme, can be sample-specific or general
#' @param u cutoff for upper extreme, can be sample-specific or general
#' @return Maximum likelihood estimates of the model parameters,
#' with 95 percent confidence intervals
#' \item{coefficients}{a vector of maximum likelihood estimates of
#' the coefficients}
#' \item{ci}{a matrix of confidence intervals for the coefficients}
#' \item{sigma}{the estimated standard deviation}
#' @details
#' The \code{\link[stats]{formula}} object is of the type
#' y~xe+xg, which describes a regression model, y=a+be*xe+bg*xg+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The covariate xg is a SNP (single-nucleotide polymorphism).
#' The variables are taken from the environment that the
#' function is called from.
#' Both xe and xg can be matrices.
#'
#' The EPS-only design is such that data is only available
#' for individuals with high and low values of the phenotype \code{y}. The
#' cut-offs \code{l} and \code{u} that specify the sampling must be given
#'
#' Thus, the data set must consist of observations only for extreme individuals
#' \code{y > u} or \code{y < l} where \code{y} is the response (phenotype)
#' of the linear regression model.
#'
#' @import MASS stats
#' @export
#' @examples
#'
#'

epsCC.lm = function(formula,l,u){
    if(class(formula)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata = model.frame(formula)
        covariates = as.matrix(model.matrix(formula)[,-1])
    options(na.action="na.omit")
    if(sum(is.na(epsdata)) > 1){stop("NA values in the data not allowed")}


    if(dim(covariates)[2]>0){
        if(dim(covariates)[2]==1){
            covnames = colnames(epsdata)[2]
        }else{
            covnames = colnames(covariates)
        }

        n = dim(epsdata)[1]
        y = epsdata[,1]

        modeldata = cbind(epsdata[,1],covariates)

        model = epsCC.loglikmax(modeldata,l,u,hessian = TRUE)

        hessian = model[[1]]
        info = -1*ginv(hessian)
        params = model[[2]]
        sigma = params[length(params)]
        coef = c()
        ci = data.frame(matrix(NA,nrow = (length(params)-1), ncol = 2))
        for(i in 1:(length(params)-1)){
            coef[i] = params[i]
            ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
            ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
        }
        colnames(ci) = c("lower 95% ci", "upper 95% ci")
        names(coef) = c("(intercept)",covnames)
        rownames(ci) = c("(intercept)",covnames)

        result = list(coef,ci,sigma)
        names(result) = c("coefficients","ci","sigma")
        return(result)
    }else{
        n = dim(epsdata)[1]
        y = epsdata[,1]
        modeldata = as.matrix(y)

        model = epsCC.loglikmax(modeldata,l,u,hessian = TRUE)

        hessian = model[[1]]
        info = -1*ginv(hessian)
        params = model[[2]]
        sigma = params[length(params)]
        coef = c()
        ci = data.frame(matrix(NA,nrow = (length(params)-1), ncol = 2))
        for(i in 1:(length(params)-1)){
            coef[i] = params[i]
            ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
            ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
        }
        colnames(ci) = c("lower 95% ci", "upper 95% ci")
        names(coef) = c("(intercept)")
        rownames(ci) = c("(intercept)")

        result = list(coef,ci,sigma)
        names(result) = c("coefficients","ci","sigma")
        return(result)
    }


}
