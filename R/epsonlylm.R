#' Fit linear model to EPS-only data
#' @description
#' \code{epsonly.lm} fits a normal linear regression model to EPS-only data
#' @param formula an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model to be fitted
#' @param cutoffs a vector \code{c(l,u)} of the lower and upper cut-offs used
#' for extreme sampling
#' @return Maximum likelihood estimates of the model parameters,
#' with 95 percent confidence intervals
#' \item{coefficients}{a vector of maximum likelihood estimates of
#' the coefficients}
#' \item{ci}{a matrix of confidence intervals for the coefficients}
#' \item{sigma}{the estimated standard deviation}
#' @details
#' The formula object is similar to that of the \code{lm} function,
#' and describes a regression model, assuming a normal distribution for
#' the residuals. See Examples.
#'
#' The data set must consist of observations only for extreme individuals
#' \code{y > u} or \code{y < l} where \code{y} is the response (phenotype)
#' of the linear regression model.
#'
#' @import MASS stats
#' @export
#' @examples
#' ## 1st example:
#' ## Create dataset:
#' N = 2000
#' xe = rnorm(n = N, mean = 2, sd = 1)
#' maf = 0.2
#' xg = sample(c(0,1,2),N,c((1-maf)^2,2*maf*(1-maf),maf^2), replace = TRUE)
#' a = 50; be = 5; bg = 0.3; sigma = 2
#' y = rnorm(N, mean = a + be*xe + bg*xg, sd = sigma)
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' y = y[extreme]
#' xe = xe[extreme]
#' xg = xg[extreme]
#' ## Fit linear model
#' epsonly.lm(y~xe+xg,cutoffs = c(l,u))

epsonly.lm = function(formula,cutoffs){
    if(length(cutoffs) != 2){stop("Invalid cutoffs vector given")}
    options(na.action="na.pass")
    epsdata = model.frame(formula)
    covariates = model.matrix(formula)[,-1]
    options(na.action="na.omit")

    n = dim(epsdata)[1]
    y = epsdata[,1]

    modelnames = attr(terms(formula), "term.labels")
    covariateorder = attr(terms(formula), "order")

    if(length(modelnames)!= dim(covariates)[2]){
        toformula = c()
        for(i in 1:length(modelnames)){
            if(covariateorder[i] > 1){
                toformula[length(toformula)+1] = modelnames[i]
            }else{
                mat = as.matrix(get(all.vars(formula)[(i+1)]))
                if(dim(mat)[2]>1){
                    for(j in 1:dim(mat)[2]){
                        assign(colnames(mat)[j],mat[,j])
                        toformula[length(toformula)+1] = colnames(mat)[j]
                    }
                }else{
                    toformula[length(toformula)+1] = modelnames[i]
                }
            }

        }
        # then there is a covariate in the fomula that is a matrix
        toformula = unique(toformula)
        formula = as.formula(paste("y ~ ", paste(toformula, collapse= "+")))
        options(na.action="na.pass")
        epsdata = model.frame(formula)
        covariates = model.matrix(formula)[,-1]
        options(na.action="na.omit")
        modelnames = attr(terms(formula), "term.labels")
    }

    if(sum(is.na(epsdata)) > 1){stop("NA values in the data not allowed")}

    modeldata = cbind(epsdata[,1],covariates)

    # Common variables in all methods, make covariates into matrices
    l = min(cutoffs)
    u = max(cutoffs)

    model = epsonlyloglikmax(modeldata,cutoffs,hessian = TRUE)

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
    names(coef) = c("(intercept)",colnames(covariates))
    rownames(ci) = c("(intercept)",colnames(covariates))

    result = list(coef,ci,sigma)
    names(result) = c("coefficients","ci","sigma")
    return(result)
}
