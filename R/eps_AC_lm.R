#' Fit linear model to EPS-full data
#' @description Fit linear model to EPS-full data
#' @param formula an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model to be fitted, see details
#' @param confounder confounding (non-genetic) covariates
#' @return Maximum likelihood estimates of the model parameters,
#' with 95 percent confidence intervals
#' \item{coefficients}{a vector of maximum likelihood estimates of the coefficients}
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
#' The EPS all case design is such that the SNP genotype is only observed
#' for individuals with high and low values of the phenotype \code{y}.
#' For remaining individuals, the unobserved genotype most be coded as NA.
#'
#' Confounders are discrete covariates (xe) and the distribution of xg is
#' modelled for each level of unique value of xe.
#'
#' @import MASS stats
#' @export
#' @examples
#'

epsAC.lm = function(formula,confounder){

    if(class(formula)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
    epsdata = model.frame(formula)
    covariates = as.matrix(model.matrix(formula)[,-1])
    options(na.action="na.omit")

    if(dim(covariates)[2]==1){covnames = colnames(epsdata)[2]
    }else{covnames = colnames(covariates)
    }

    n = dim(epsdata)[1]
    y = epsdata[,1]
    if(sum(is.na(y))>0){
        stop("NA-values not allowed in response variable of regression model")}

    # Which covariates are snps, covariates, interactions?
    # Interactions first
    interact = FALSE # Refers to gene-environ interactions
    interactind = list()
    xxid = c()
    xsnpid = c()

    modelnames = attr(terms(formula), "term.labels")
    covariateorder = attr(terms(formula), "order")

    nint = sum(covariateorder>1)
    if(nint > 0){
        covintnames = modelnames[covariateorder > 1]
        for(i in 1:nint){
            if((!strsplit(covintnames[i],":")[[1]][1]%in%modelnames) |
               (!strsplit(covintnames[i],":")[[1]][2]%in%modelnames)){
                stop("No interactions without corresponding main effect in the model")
            }
            intid = which(covnames == covintnames[i])
            if(sum(is.na(covariates[,intid]))>0){
                # Interaction with SNP
                interact = TRUE
                xsnpid[length(xsnpid)+1] = intid
                # also collect ids in list for other functions to use
                term = covintnames[i]
                t = match(strsplit(term,":")[[1]], covnames)
                interactind[[(length(interactind)+1)]] = t
            }else{
                xxid[length(xxid)+1] = intid
            }
        }
    }

    # we have located which columns in covariates are interactions
    # now we classify main environment and snp
    snpid = c()
    xid = c()
    colstocheck = c(1:dim(covariates)[2])
    if(length(xsnpid)>0 | length(xxid)>0){
        colstocheck = colstocheck[-c(xxid,xsnpid)]
    }
    for(i in 1:length(colstocheck)){
        if(sum(is.na(covariates[,colstocheck[i]]))>0){
            snpid[length(snpid)+1] = colstocheck[i]
        }else{
            xid[length(xid)+1] = colstocheck[i]
        }
    }

    if(interact){
        for(i in 1:length(interactind)){
            t = interactind[[i]]
            if(t[1] %in% snpid){
                # always genetic variant first
                interactind[[i]] = c(match(t[1],snpid),match(t[2],xid))
            }else{ interactind[[i]] = c(match(t[2],snpid),match(t[1],xid)) }
        }
    }

    message(paste("Identified ", length(snpid), " SNP(s), ",
                  length(xid), " environmental covariate(s), ",
                  length(xxid), " environmental interaction(s), and ",
                  length(xsnpid), " gene-environment interaction(s). ",
                  sep = ""))

    xg = as.matrix(covariates[,snpid])
    isxe = FALSE
    if(length(xid) > 0){
        isxe = TRUE
        xe = as.matrix(covariates[,c(xid,xxid)])
    }
    newnames = c(covnames[c(xid,xxid)],covnames[snpid],covnames[xsnpid])

    if(!interact){
        #################################################################
        # No interactions
        #################################################################
        if(!isxe){
            ###############################################################
            #  No environmental covariates
            ###############################################################
            data = cbind(y,xg)
            ng = dim(as.matrix(xg))[2]
            model = epsAC.loglikmax(data,ng,hessian = TRUE)
            params = model[[2]]
            nparam = 1 + ng
            sigma = params[(nparam + 1)]
            hessian = model[[1]]
            info = -1*ginv(hessian)
            coef = c()
            ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
            for(i in 1:nparam){
                coef[i] = params[i]
                ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
            }
            gfreqs = model[[3]]
            genotypes = model[[4]]
            resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
            colnames(resg) = c("P(Xg)",covnames[snpid])
            colnames(ci) = c("lower 95% ci", "upper 95% ci")
            names(coef) = c("(intercept)",newnames)
            rownames(ci) = c("(intercept)",newnames)
            result = list(coef,ci,sigma,resg)
            names(result) = c("coefficients","ci","sigma","Xg")
            return(result)
        }else{
            ###############################################################
            # Environmental covariates
            ###############################################################
            data = cbind(y,xe,xg)
            ng = dim(as.matrix(xg))[2]
            ne = dim(as.matrix(xe))[2]

            if(missing(confounder)){
                ###########################################################
                # No confounders
                ###########################################################
                model = epsAC.loglikmax(data,ng,hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ne + ng # alpha, betae, betag
                sigma = params[(nparam + 1)]
                hessian = model[[1]]
                info = -1*ginv(hessian)
                coef = c()
                ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
                for(i in 1:(1 + ne + ng)){
                    coef[i] = params[i]
                    ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                    ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                }
                gfreqs = model[[3]]
                genotypes = model[[4]]
                resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
                colnames(resg) = c("P(Xg)",covnames[snpid])
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",newnames)
                rownames(ci) = c("(intercept)",newnames)
                result = list(coef,ci,sigma,resg)
                names(result) = c("coefficients","ci","sigma","Xg")
                return(result)
            }else{
                #########################################################
                # Confounders
                #########################################################
                model = epsAC.loglikmaxcond(data, confounder, ng,
                                            hessian = TRUE,snpnames=covnames[snpid])
                params = model[[2]]
                nparam = 1 + ne + ng
                sigma = params[(nparam + 1)]
                probs = params[(nparam + 2):length(params)]
                hessian = model[[1]]
                info = -1*ginv(hessian)
                coef = c()
                ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
                for(i in 1:(1 + ne + ng)){
                    coef[i] = params[i]
                    ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                    ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                }
                resg = model[[3]]
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",
                                newnames)
                rownames(ci) = c("(intercept)",
                                 newnames)
                result = list(coef,ci,sigma,resg)
                names(result) = c("coefficients","ci","sigma","Xg")
                return(result)
            }
        }
    }else{
        #################################################################
        #  Interactions
        #################################################################

        # stop("Interaction analysis currently not available")
        #
        # data = cbind(y,xe,xg)
        # xe = as.matrix(xe)
        # xg = as.matrix(xg)
        # ng = dim(as.matrix(xg))[2]
        # ne = dim(as.matrix(xe))[2]
        # intnames = c()
        # for(i in 1:length(interactind)){
        #     intnames[i] = paste(colnames(xe)[interactind[[i]][2]],"*",
        #                         colnames(xg)[interactind[[i]][1]],sep="")
        # }
        # if(!confounder){
        #     ###############################################################
        #     # 2.1.1: No confounders
        #     ###############################################################
        #     model = epsAC.loglikmaxint(data, ng, interactind = interactind,hessian = TRUE)
        #     params = model[[2]]
        #     nparam = 1 + ne + ng + length(interactind)
        #     sigma = params[(nparam + 1)]
        #     hessian = model[[1]]
        #     info = -1*ginv(hessian)
        #     coef = c()
        #     ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
        #     for(i in 1:nparam){
        #         coef[i] = params[i]
        #         ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
        #         ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
        #     }
        #     gfreqs = model[[3]]
        #     genotypes = model[[4]]
        #     resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
        #     colnames(resg) = c("P(Xg)",covnames[snpid])
        #     colnames(ci) = c("lower 95% ci", "upper 95% ci")
        #     names(coef) = c("(intercept)",newnames)
        #     rownames(ci) = c("(intercept)",newnames)
        #     result = list(coef,ci,sigma,resg)
        #     names(result) = c("coefficients","ci","sigma","Xg")
        #     return(result)
        # }else{
        #     ###############################################################
        #     # 2.1.1: Confounders
        #     ###############################################################
        #     model = epsAC.loglikmaxcondint(data, ng, cind = cind,
        #                                    interactind = interactind,
        #                                    hessian = TRUE,snpnames=covnames[snpid])
        #     params = model[[2]]
        #     nparam = 1 + ne + ng + length(interactind)
        #     sigma = params[(nparam + 1)]
        #     hessian = model[[1]]
        #     info = -1*ginv(hessian)
        #     coef = c()
        #     ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
        #     for(i in 1:nparam){
        #         coef[i] = params[i]
        #         ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
        #         ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
        #     }
        #     resg = model[[3]]
        #     colnames(ci) = c("lower 95% ci", "upper 95% ci")
        #     names(coef) = c("(intercept)",
        #                     newnames)
        #     rownames(ci) = c("(intercept)",
        #                      newnames)
        #     result = list(coef,ci,sigma,resg)
        #     names(result) = c("coefficients","ci","sigma","Xg")
        #     return(result)
        # }
    }
}
