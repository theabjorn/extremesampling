#' @title Test for associations under the EPS complete-case design
#' @description
#' \code{epsAC.rv.test} performs a score test for common genetic variants
#' under the EPS complete-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param RV a matrix of genetic variants to be tested against the null
#' @return \code{epsAC.rv.test} returns for each genetic variant:
#' \item{statistic}{the score test statistic}
#' \item{p.value}{the P-value}
#' @details
#' nanana
#'
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

epsAC.rv.test = function(nullmodel,RV,confounder = FALSE, cx,method = "naive",weights){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    isxe = FALSE
    if(dim(covariates0)[2]>0){
        isxe = TRUE
        xe = covariates0
        if(dim(covariates0)[2] ==1){
            colnames(covariates0) = colnames(epsdata0)[2]
        }
        # Confounders
        if(confounder){
            if(missing(cx)){cx = colnames(covariates0)
            }else if(is.na(match(cx,colnames(covariates0)))){
                stop("The name of the confounder given is not the name of a covariate.")
            }
            message(paste("Confounders: ", toString(cx),sep=""))
            cind = match(cx,colnames(covariates0))
            xecind = as.matrix(covariates0[,cind])
            for(j in 1:length(cx)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                         Please recode you confounder to satisfy this.")
                }
                }
            }
        }

    RV = as.matrix(RV)
    #if(dim(RV)[2]<2){
    #    stop("RV must be a set of more than one genetic variant")
    #}

    if(missing(weights)){
        weights = rep(1,dim(RV)[2])
    }

    if(method == "naive"){
        epsAC.rv.test.naive(epsdata0,covariates0,RV,isxe,confounder,cind)
    }else if(method == "burden"){
        g = rep(0,dim(RV)[1])
        for(c in 1:dim(RV)[2]){
            g = g + c(weights[c])*c(RV[,c])
        }
        epsAC.test(nullmodel,SNP=g,confounder,cx)
    }else if(method == "lmm"){
        epsAC.rv.test.lmm(epsdata0,covariates0,RV,isxe,confounder,cind,weights)
    }
}

