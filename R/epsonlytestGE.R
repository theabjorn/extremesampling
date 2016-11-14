#' @title Test gene-environment interactions in the EPS-only design
#' @description
#' \code{epsonly.testGE} performs a score test for gene-environment interactions
#' in the EPS-only design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param GE a list of interactions, with the colon-symbol used to denote
#' interaction, e.g. \code{c("xe1:xg1","xe1:xg2")}
#' @param cutoffs a vector \code{c(l,u)} of the lower and upper cut-offs used
#' for extreme sampling
#' @param onebyone \code{TRUE} if interactions should be tested one by one
#' for inclusion in the model, default \code{TRUE}
#' @return \code{epsonly.testGE} returns
#' \item{statistic}{the value of the score test statistic}
#' \item{parameter}{the degrees of freedom of the statistic}
#' \item{p.value}{the P-value for the test}
#' @details
#' The \code{nullmodel} formula object are similar to that of the
#' \code{lm} function, and describes a regression model,
#' assuming a normal distribution for the residuals. The data are taken
#' from the environment from which the function is called.
#'
#' If there are more than one additional gene-environment interaction
#' in \code{GE} to be tested against the null model, then if one-by-one is TRUE,
#' each additional variable is seperately tested against the null model.
#'
#' The data set must consist of observations only for extreme individuals
#' \code{y > u} or \code{y < l} where \code{y} is the response (phenotype)
#' of the linear regression model.
#'
#' @import MASS stats
#' @export
#' @examples
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
#' ## Perform score test:
#' epsonly.testGE(y~xe+xg,c("xe:xg"),cutoffs = c(l,u))

epsonly.testGE = function(nullmodel,GE,cutoffs,onebyone = TRUE){
    if(length(cutoffs) != 2){stop("Invalid cutoffs vector given")}

    options(na.action="na.pass")
    epsdata0 = model.frame(nullmodel)
    covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    if(sum(is.na(epsdata0) > 1)){
        stop("NA values in the data not allowed")}

    modeldata = cbind(epsdata0[,1],covariates0)
    l = min(cutoffs)
    u = max(cutoffs)

    covariateorder = attr(terms(nullmodel), "order")
    maineffects = attr(terms(nullmodel),"term.labels")[covariateorder == 1]

    nge = length(GE)
    n = dim(epsdata0)[1]
    xge = matrix(NA,ncol = nge,nrow = n)
    for(i in 1:nge){
        t = match(strsplit(GE[i],":")[[1]], maineffects)
        if(sum(is.na(t)>0)){
            stop(paste("Cannot include interaction term ",toString(GE[i]),
                       ", because corresponding main effects were not found in the null model",
                       sep = ""))
        }
        xge[,i] = covariates0[,t[1]]*covariates0[,t[2]]
    }

    y = epsdata0[,1]

    if(onebyone){
        ###################################################################
        # Test gene-environment interactions one by one
        ###################################################################

        x = as.matrix(covariates0)

        fit = epsonlyloglikmax(modeldata,c(l,u)) # Fit under H0
        alpha = fit[1]
        beta = fit[2:(length(fit)-1)]
        sigma = fit[length(fit)]
        sigma2 = sigma*sigma
        xbeta = x%*%beta
        zl = (l-alpha-xbeta)/sigma
        zu = (u-alpha-xbeta)/sigma

        h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
        h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
        h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
        h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
            (1-pnorm(zu)+pnorm(zl))

        a = 1 - h1 - h0*h0
        b = h0 - h2 - h0*h1
        c = 2 + 2*h1 - h3 - h1*h1

        I11_11 = sum(a)
        I11_22 = t(x)%*%(diag(a[,1])%*%x)

        I11_21 = t(t(a)%*%x)
        I11_12 = t(a)%*%x

        I11_33 = sum(c)
        I11_31 = sum(b)
        I11_13 = sum(b)
        I11_23 = t(t(b)%*%x)
        I11_32 = t(b)%*%x

        I11 = cbind(rbind(I11_11,I11_21,I11_31),
                    rbind(I11_12,I11_22,I11_32),
                    rbind(I11_13,I11_23,I11_33))

        statistic = matrix(NA,ncol = 1, nrow = nge)
        parameter = matrix(NA,ncol = 1, nrow = nge)
        pvalue = matrix(NA,ncol = 1, nrow = nge)
        rownames(statistic) = GE
        rownames(parameter) = GE
        rownames(pvalue) = GE

        colnames(statistic) = "t"
        colnames(parameter) = "d.o.f"
        colnames(pvalue) = "p.value"

        for(i in 1:nge){
            gi = xge[,i]
            I22 = t(gi)%*%(diag(a[,1])%*%gi)

            I21_1 = t(t(a)%*%gi)
            I21_2 = t(t(x)%*%(diag(a[,1])%*%gi))
            I21_3 = t(b)%*%gi
            I21 = cbind(I21_1,I21_2,I21_3)
            I12 = t(I21)

            Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
            s = sum((y-alpha - xbeta +sigma*h0)*gi)/sigma2
            t = s*s/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            statistic[i,] = t
            parameter[i,] = 1
            pvalue[i,] = pval
        }
        result = list(statistic,parameter,pvalue)
        names(result) = c("statistic","parameter","p.value")
        return(result)

    }else{
        ###################################################################
        # Test all interactions a time
        ###################################################################
        if(nge > 10){
            stop("Do not test more than 10 interactions simultaneously,
                 choose onebyone = TRUE")
        }
        x = as.matrix(covariates0)

        fit = epsonlyloglikmax(modeldata,c(l,u))
        alpha = fit[1]
        beta = fit[2:(length(fit)-1)]
        sigma = fit[length(fit)]
        sigma2 = sigma*sigma

        xbeta = x%*%beta
        zl = (l-alpha-xbeta)/sigma
        zu = (u-alpha-xbeta)/sigma

        h0 = (-dnorm(zu)+dnorm(zl))/(1-pnorm(zu)+pnorm(zl))
        h1 = (-dnorm(zu)*zu+dnorm(zl)*zl)/(1-pnorm(zu)+pnorm(zl))
        h2 = (-dnorm(zu)*zu*zu+dnorm(zl)*zl*zl)/(1-pnorm(zu)+pnorm(zl))
        h3 = (-dnorm(zu)*zu*zu*zu+dnorm(zl)*zl*zl*zl)/
            (1-pnorm(zu)+pnorm(zl))

        a = 1 - h1 - h0*h0
        b = h0 - h2 - h0*h1
        c = 2 + 2*h1 - h3 - h1*h1

        I11_11 = sum(a)
        I11_22 = t(x)%*%(diag(a[,1])%*%x)

        I11_21 = t(t(a)%*%x)
        I11_12 = t(a)%*%x

        I11_33 = sum(c)
        I11_31 = sum(b)
        I11_13 = sum(b)
        I11_23 = t(t(b)%*%x)
        I11_32 = t(b)%*%x

        I11 = cbind(rbind(I11_11,I11_21,I11_31),
                    rbind(I11_12,I11_22,I11_32),
                    rbind(I11_13,I11_23,I11_33))

        I22 = t(xge)%*%(diag(a[,1])%*%xge)

        I21_1 = t(t(a)%*%xge)
        I21_2 = t(t(x)%*%(diag(a[,1])%*%xge))
        I21_3 = t(t(b)%*%xge)
        I21 = cbind(I21_1,I21_2,I21_3)
        I12 = t(I21)

        Sigma = (1/sigma2)*(I22 - I21%*%ginv(I11)%*%t(I21))
        s = t(y-alpha - xbeta +sigma*h0)%*%xge/sigma2
        t = s%*%ginv(Sigma)%*%t(s)
        pval = pchisq(t,nge,lower.tail=FALSE)
        result = list(c(t),nge,c(pval))
        names(result) = c("statistic","parameter","p.value")
        return(result)
        }
}
