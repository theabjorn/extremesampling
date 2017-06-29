#' @title Score test EPS-CC
#' @description
#' \code{epsCC.test} performs a score test for common genetic variants
#' under the EPS all-case design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis, see details
#' @param SNP a matrix of variables to be tested against the null (NA for
#' not genotyped individuals)
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @return \code{epsfull.test} returns
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
#' The covariate \code{xg} is a SNP (single-nucleotide polymorphism). If
#' SNP is a matrix, each SNP (column) is tested against the null model.
#'
#' For the EPS all-case design, the SNP genotypes are only observed
#' for individuals with extreme values of the phenotype \code{y}, and possibly
#' some random samples. For remaining individuals, the unobserved genotype
#' must be coded as NA.
#'
#' If confounder = TRUE, the genetic variables are assumed to have
#' different distribution for different levels of other (non-genetic)
#' covariates, these can be specified by a vector of names \code{cx}.
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
#' extreme = (y < l) | (y >= u)
#' # Create the EPS-full data set by setting
#' # the SNP values of non-extremes to NA
#' xg1[!extreme] = NA; xg2[!extreme] = NA; xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#' # Testing
#' epsAC.test(y~xe1+xe2,SNP = xg)$p.value

epsAC.test = function(nullmodel, SNP, confounder = FALSE, cx){
    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    if(is.null(colnames(SNP))){
        SNP = as.matrix(SNP)
        colnames(SNP) = paste0("SNP",1:dim(SNP)[2])}

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

    totest = colnames(SNP)
    xg = as.matrix(SNP)
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

    if(isxe & !confounder){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # No confounding assumed
        ###############################################################
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        coefs = fit$coef
        alpha = coefs[1]
        beta = coefs[2:length(coefs)]
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

            egvec = rep(eg,n)

            sig1 = (t(egvec)%*%(Ie-He)%*%egvec)/sigma2
            sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

            Sigma = sig1 + sig2

            s = (sum(z[extreme]*g) +sum(z[!extreme]*eg))/sigma2
            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            statistic[i,] = t
            pvalue[i,] = pval
        }
        result = list(statistic,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }else if(isxe & confounder){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # Confounding assumed
        ###############################################################
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        coefs = fit$coef
        alpha = coefs[1]
        beta = coefs[2:length(coefs)]
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
            uindex = list()
            for(u in 1:nu){
                uindex[[u]] = (xe[extreme,cind] == ux[u,])
            }

            eg = c()
            varg = c()
            egvec = rep(0,n)

            for(u in 1:nu){
                uind = uindex[[u]]
                allind = (xe[,cind] == ux[u,])

                egvec[allind] = mean(g[uind])

                eg[u] = mean(g[uind])
                varg[u] = var(g[uind])
            }

            sig1 = (t(egvec)%*%(Ie-He)%*%egvec)/sigma2
            sig2 = 0

            s = sum(z[extreme]*g)

            for(u in 1:nu){
                uind = uindex[[u]]
                n_uc = length(z[extreme][uind])
                sig2 = sig2 + varg[u]*(sum(z[extreme][uind]^2) - (1/n_uc)*(sum(z[extreme][uind]))^2)/(sigma2^2)

                tempind = (xe[!extreme,cind] == ux[u,])
                s = s + sum(z[!extreme][tempind]*eg[u])
            }

            s = s/sigma2
            Sigma = sig1 + sig2

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

        for(i in 1:ng){
            extreme = !is.na(xg[,i])
            g = as.matrix(xg[extreme,i])
            n_cc = sum(extreme)

            eg = mean(g)
            varg = var(g)

            sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

            Sigma = sig2

            s = (sum(z[extreme]*g) +sum(z[!extreme]*eg))/sigma2
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
