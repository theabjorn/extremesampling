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
#' N = 5000 # Number of individuals in a population
#' xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
#' xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
#' xg1 = sample(c(0,1,2),N,c(0.4,0.3,0.3), replace = TRUE) # xg
#' xg2 = sample(c(0,1,2),N,c(0.5,0.3,0.2), replace = TRUE) # xg
#' # Model parameters
#' a = 50; be1 = 5; be2 = 8; bg1 = 0.3; bg2 = 0.6; sigma = 2
#' # Generate response y
#' y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 + bg2*xg2, sd = sigma)
#' # Identify extremes, here upper and lower 25% of population
#' u = quantile(y,probs = 3/4,na.rm=TRUE)
#' l = quantile(y,probs = 1/4,na.rm=TRUE)
#' extreme = (y < l) | (y >= u)
#' # Create the EPS-full data set by setting
#' # the xg values of non-extremes to NA
#' xg1[!extreme] = NA; xg2[!extreme] = NA; xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#' # Testing
#' epsAC.test(y~xe1+xe2,xg = xg)$p.value

epsAC.test = function(nullmodel, xg, confounder){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(is.null(colnames(xg))){
        xg = as.matrix(xg)
        colnames(xg) = paste0("xg",1:dim(xg)[2])}

    options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    conf = FALSE
    if(!missing(confounder)){
        conf = TRUE
        xec = as.matrix(confounder)
        message(paste("Confounding assumed"))
        for(j in 1:dim(xec)[2]){
        if(length(unique(xec[,j])) > 10){
            stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                 Please recode you confounder to satisfy this.")
            }
        }
    }

    isxe = FALSE
    if(dim(covariates0)[2]>0){
        isxe = TRUE
        xe = covariates0
        if(dim(covariates0)[2] ==1){
            colnames(covariates0) = colnames(epsdata0)[2]
        }
    }

    totest = colnames(xg)
    xg = as.matrix(xg)
    ng = dim(xg)[2]


    ###################################################################
    # Test genetic covariates one by one
    ###################################################################
    statistic = matrix(NA,ncol = 1, nrow = ng)
    pvalue = matrix(NA,ncol = 1, nrow = ng)
    rownames(statistic) = totest
    rownames(pvalue) = totest
    colnames(statistic) = "t"
    colnames(pvalue) = "p.value"

    if(isxe & !conf){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # No confounding assumed
        ###############################################################
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)
        Xe = cbind(rep(1,n),xe)

        meanimpute = function(g){
            g[is.na(g)] = mean(g,na.rm=TRUE)
            return(g)
        }
        xgm = apply(xg,2,meanimpute)
        veg = crossprod(Xe,xgm)
        veeinvveg = solve(crossprod(Xe,Xe))%*%veg
        var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2
        calcvar2 = function(g){
            varg = var(g,na.rm=TRUE)
            n_ic = sum(is.na(g))
            n_cc = sum(!is.na(g))
            return(varg*(sum(z[is.na(g)]^2) - sigma2*n_ic + (1/n_cc)*(sum(z[is.na(g)]))^2)/(sigma2^2))
        }
        var2 = apply(xg,2,calcvar2)

        s = c(t(z)%*%xgm/sigma2)
        Sigma = var1 - var2
        t = (s*s)/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        statistic = matrix(t,ncol = 1, nrow = ng)
        pvalue = matrix(pval,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)

    }else if(isxe & conf){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # Confounding assumed
        ###############################################################
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = cbind(rep(1,n),xe)

        ux = as.matrix(unique(xec))
        nu = dim(ux)[1]

        meanimpute = function(g){
            for(u in 1:nu){
                uindex_all = (xec == ux[u,])
                gu = g[uindex_all]
                gu[is.na(gu)] = mean(gu,na.rm=TRUE)
                g[uindex_all] = gu
            }
            return(g)
        }

        xgm = apply(xg,2,meanimpute)
        veg = crossprod(Xe,xgm)
        veeinvveg = solve(crossprod(Xe,Xe))%*%veg
        var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2

        calcvar2 = function(g){
            tmp = 0
            for(u in 1:nu){
                uindex_all = (xec == ux[u,])
                gu = g[uindex_all]
                varg = var(gu,na.rm=TRUE)
                zu = z[uindex_all]
                n_uic = sum(is.na(gu))
                n_ucc = sum(!is.na(gu))
                tmp = tmp + varg*(sum(zu[is.na(gu)]^2) - sigma2*n_uic + (1/n_ucc)*(sum(zu[is.na(gu)]))^2)/(sigma2^2)
            }
            return(tmp)
        }
        var2 = apply(xg,2,calcvar2)

        s = c(t(z)%*%xgm/sigma2)
        Sigma = var1 - var2
        t = (s*s)/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        statistic = matrix(t,ncol = 1, nrow = ng)
        pvalue = matrix(pval,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }else if(conf){
        fit = lm(y~1) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = matrix(1,ncol=1,nrow = n)

        ux = as.matrix(unique(xec))
        nu = dim(ux)[1]

        meanimpute = function(g){
            for(u in 1:nu){
                uindex_all = (xec == ux[u,])
                gu = g[uindex_all]
                gu[is.na(gu)] = mean(gu,na.rm=TRUE)
                g[uindex_all] = gu
            }
            return(g)
        }

        xgm = apply(xg,2,meanimpute)
        veg = crossprod(Xe,xgm)
        veeinvveg = solve(crossprod(Xe,Xe))%*%veg
        var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2

        calcvar2 = function(g){
            tmp = 0
            for(u in 1:nu){
                uindex_all = (xec == ux[u,])
                gu = g[uindex_all]
                varg = var(gu,na.rm=TRUE)
                zu = z[uindex_all]
                n_uic = sum(is.na(gu))
                n_ucc = sum(!is.na(gu))
                tmp = tmp + varg*(sum(zu[is.na(gu)]^2) - sigma2*n_uic + (1/n_ucc)*(sum(zu[is.na(gu)]))^2)/(sigma2^2)
            }
            return(tmp)
        }
        var2 = apply(xg,2,calcvar2)

        #s = c(t(z[!is.na(xg[,1])])%*%xg[!is.na(xg[,1]),]/sigma2)
        s = c(t(z)%*%xgm/sigma2)
        Sigma = var1 - var2
        t = (s*s)/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        statistic = matrix(t,ncol = 1, nrow = ng)
        pvalue = matrix(pval,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }else{
        ##############################################################
        # no environmental covariates (xe) present in the null model
        ##############################################################
        fit = lm(y~1) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = matrix(1,ncol=1,nrow = n)

        meanimpute = function(g){
            g[is.na(g)] = mean(g,na.rm=TRUE)
            return(g)
        }
        xgm = apply(xg,2,meanimpute)
        veg = crossprod(Xe,xgm)
        veeinvveg = solve(crossprod(Xe,Xe))%*%veg
        var1 = (colSums(xgm*xgm)-colSums(veg*veeinvveg))/sigma2
        calcvar2 = function(g){
            varg = var(g,na.rm=TRUE)
            n_ic = sum(is.na(g))
            n_cc = sum(!is.na(g))
            return(varg*(sum(z[is.na(g)]^2) - sigma2*n_ic + (1/n_cc)*(sum(z[is.na(g)]))^2)/(sigma2^2))
        }
        var2 = apply(xg,2,calcvar2)

        s = c(t(z)%*%xgm/sigma2)
        #s = c(t(z[!is.na(xg[,1])])%*%xg[!is.na(xg[,1]),]/sigma2)
        Sigma = var1 - var2
        t = (s*s)/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        statistic = matrix(t,ncol = 1, nrow = ng)
        pvalue = matrix(pval,ncol = 1, nrow = ng)
        rownames(statistic) = totest
        rownames(pvalue) = totest
        colnames(statistic) = "t"
        colnames(pvalue) = "p.value"
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)

    }

}
