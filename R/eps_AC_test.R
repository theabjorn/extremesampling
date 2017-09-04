#' @title Score test EPS-AC
#' @description
#' \code{epsAC.test} performs a score test for common genetic variants
#' under the EPS all-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis, see details
#' @param xg a matrix of variables to be tested against the null (NA for
#' not genotyped individuals)
#' @param confounder (optional) vector of names of confounding (non-genetic) covariates
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
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        for(i in 1:ng){
            extreme = !is.na(xg[,i])
            g = as.matrix(xg[extreme,i])
            n_cc = sum(extreme)

            eg = mean(g)
            varg = var(g)

            gm = rep(eg,n)
            gm[extreme] = g

            var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
            var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

            Sigma = var1 - var2

            s = (t(gm)%*%(Ie-He)%*%y)/sigma2

            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            statistic[i,] = t
            pvalue[i,] = pval
        }
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
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        for(i in 1:ng){
            extreme = !is.na(xg[,i])
            g = as.matrix(xg[extreme,i])
            n_cc = sum(extreme)

            # confounding
            ux = as.matrix(unique(xe[extreme,cind]))
            nu = dim(ux)[1]
            if(nu != dim(as.matrix(unique(xe[,cind])))[1]){
                warning("All unique levels of confounder not found in extreme sample")
            }
            uindex_all = list()
            uindex_cc = list()
            uindex_ic = list()
            for(u in 1:nu){
                uindex_all[[u]] = (xe[,cind] == ux[u,])
                uindex_cc[[u]] = (xe[extreme,cind] == ux[u,])
                uindex_ic[[u]] = (xe[!extreme,cind] == ux[u,])
            }

            eg = c()
            varg = c()
            egvec = rep(0,n)

            for(u in 1:nu){
                uind = uindex_cc[[u]]
                allind = uindex_all[[u]]
                egvec[allind] = mean(g[uind])
                eg[u] = mean(g[uind])
                varg[u] = var(g[uind])
            }

            gm = rep(0,n)
            gm[extreme] = g
            gm[!extreme] = egvec[!extreme]

            var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2

            var2 = 0
            for(u in 1:nu){
                uind_ic = uindex_ic[[u]]
                uind_cc = uindex_cc[[u]]
                uind_all = uindex_all[[u]]

                n_uic = length(z[!extreme][uind_ic])
                n_ucc = length(z[extreme][uind_cc])

                var2 = var2 + varg[u]*(sum(z[!extreme][uind_ic]^2) - sigma2*(n_uic) - (1/n_ucc)*(sum(z[!extreme][uind_ic]))^2)/(sigma2^2)
            }

            Sigma = var1 - var2
            s = (t(gm)%*%(Ie-He)%*%y)/sigma2

            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            statistic[i,] = t
            pvalue[i,] = pval
        }
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }else{
        ##############################################################
        # no environmental covariates (xe) present in the null model
        ##############################################################
        fit = lm(y~1) # Fit under H0
        z = fit$residuals
        alpha = mean(y)
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = rep(1,n)
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        for(i in 1:ng){
            extreme = !is.na(xg[,i])
            g = as.matrix(xg[extreme,i])
            n_cc = sum(extreme)

            eg = mean(g)
            varg = var(g)

            gm = rep(eg,n)
            gm[extreme] = g

            var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
            var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

            Sigma = var1 - var2

            s = (t(gm)%*%(Ie-He)%*%y)/sigma2

            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            statistic[i,] = t
            pvalue[i,] = pval
        }
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }

}
