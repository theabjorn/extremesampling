#' Fit linear model to EPS-complete data
#' @description
#' \code{epscomp.lm} fits a normal linear regression model to EPS-complete data
#' @param formula an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model to be fitted, see details
#' @param hwe \code{TRUE} if Hardy-Weinberg equilibrium is assumed, default
#' set to \code{FALSE}
#' @param maf optional value for the minor allele frequencies under HWE
#' @param gfreq frequencies of genotypes for each SNP if known or
#' otherwise estimated
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @return Maximum likelihood estimates of the model parameters,
#' with 95 percent confidence intervals
#' \item{coefficients}{a vector of maximum likelihood estimates of the coefficients}
#' \item{ci}{a matrix of confidence intervals for the coefficients}
#' \item{sigma}{the estimated standard deviation}
#' \item{maf}{the estimated minor allele frequencies, returned if
#' \code{HWE = TRUE} but \code{MAF} is unspecified}
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
#' The EPS-complete design is such that the SNP genotype is only observed
#' for individuals with high and low values of the phenotype \code{y}.
#' For remaining individuals, the unobserved genotype most be coded as NA.
#' A SNP is assumed to have possible genotype 0, 1 or 2 according to the
#' number of minor-alleles. The distribution of the genotype is assumed
#' unknown and multinomially distributed. I.e. P(xg=0) = p0, P(xg=1) = p1,
#' and P(xg=2) = p2 = 1-p0-p1.
#' Hardy-Weinberg equilibrium with known or uknown MAF can be assumed,
#' then p0 = (1-q)^2, p1 = 2q(1-q) and p2 = q^2, where q is the MAF.
#' The parameters p0, p1, p2 can be given in \code{gfreq} if they are known.
#'
#' If confounder = TRUE, the genetic variables are assumed to be
#' multinomally distributed, with different distribution
#' for different levels of other (non-genetic) covariates, these can
#' be specified by a vector of names \code{cx}.
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
#' # Create the EPS-complete data set by setting
#' # the SNP values of non-extremes to NA
#' xg1[!extreme] = NA
#' xg2[!extreme] = NA
#' xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#'
#' # Fit model
#' epscomp.lm(y~xe1+xe2+xg1+xg2)
#' # Alternatives
#' # epscomp.lm(y~xe+xg)
#' # epscomp.lm(y~xe+xg,hwe = TRUE)
#'
#' # Model with interaction term
#' epscomp.lm(y~xe+xg+xe1*xg2)
#'

epscomp.lm = function(formula, hwe = FALSE, maf, gfreq,
                      confounder = FALSE, cx){

    if(class(formula)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
        epsdata = model.frame(formula)
        covariates = model.matrix(formula)[,-1]
    options(na.action="na.omit")

    n = dim(epsdata)[1]
    y = epsdata[,1]
    if(sum(is.na(y))>0){
        stop("NA-values not allowed in response variable of regression model")}

    modelnames = attr(terms(formula), "term.labels")
    covariateorder = attr(terms(formula), "order")
    if(length(modelnames)!= dim(covariates)[2]){
        # there is a covariate in the formula that is a matrix
        toformula = c()
        # keep track of column index in covariates
        covind = 0
        for(i in 1:length(modelnames)){
            if(covariateorder[i]>1){
                covind = covind + 1
                assign(modelnames[i],covariates[,covind])
                toformula[length(toformula)+1] = modelnames[i]
            }else{
                mat = as.matrix(get(modelnames[i]))
                if(dim(mat)[2]>1){ # matrix covariate found
                    if(is.null(colnames(mat))){
                        warning("Matrix of covariates should have column names")
                        colnames(mat) = paste0(modelnames[1],1:dim(mat)[2])
                    }
                    for(j in 1:dim(mat)[2]){
                        covind = covind + 1
                        assign(colnames(mat)[j],mat[,j]) # make a new vector
                        toformula[length(toformula)+1] = colnames(mat)[j]
                    }
                }else{
                    covind = covind + 1
                    toformula[length(toformula)+1] = modelnames[i]
                }
            }
        }
        toformula = unique(toformula)
        formula = as.formula(paste("y ~ ", paste(toformula, collapse= "+")))
        options(na.action="na.pass")
        epsdata = model.frame(formula)
        covariates = model.matrix(formula)[,-1]
        options(na.action="na.omit")
        modelnames = attr(terms(formula), "term.labels")
    }

    # Which covariates are snps, covariates, interactions?
    covariateorder = attr(terms(formula), "order")
    interact = FALSE
    snpid = c()
    xid = c()
    xxid = c()
    xsnpid = c()
    for(c in 1:dim(covariates)[2]){
        if(covariateorder[c] == 1 & sum(is.na(covariates[,c])) > 0){
            snpid[length(snpid)+1] = c
        }else if(covariateorder[c] > 1 & sum(is.na(covariates[,c])) > 0){
            xsnpid[length(xsnpid)+1] = c
            interact = TRUE
        }else if(covariateorder[c] == 1){
            xid[length(xid)+1] = c
        }else if(covariateorder[c] > 1){
            xxid[length(xxid)+1] = c
        }
    }

    if(interact & is.null(snpid)){
        stop("No SNP-environment interaction without a SNP main effect")
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
    newnames = c(modelnames[c(xid,xxid)],modelnames[snpid],modelnames[xsnpid])

    if(hwe & missing(maf)){
        maf = NA
    }
    if(missing(gfreq)){gfreq = NA}

    maineffects = attr(terms(formula),"term.labels")[covariateorder == 1]
    if(interact){
        interactind = list()
        message(paste("Interactions: ",modelnames[xsnpid],sep=""))
        for(i in 1:length(xsnpid)){
            term = modelnames[xsnpid][i]
            t = match(strsplit(term,":")[[1]], maineffects)
            if(sum(is.na(t))>1){
                stop("No interactions without corresponding main effect in the model")
            }
            if(t[1] %in% snpid){
                # always genetic variant first
                interactind[[i]] = c(match(t[1],snpid),match(t[2],xid))
            }else{ interactind[[i]] = c(match(t[2],snpid),match(t[1],xid)) }
        }
    }

    if(confounder & isxe){
        if(hwe){stop("Confounders and Hardy-Weinberg equilibrium not allowed simultaneously")}
        if(missing(cx)){cx = maineffects[xid]}
        message(paste("Confounders: ", toString(cx),sep=""))
        cind = match(match(cx,maineffects),xid)
        xecind = as.matrix(covariates[,cind])
        for(j in 1:length(cx)){
            if(length(unique(xecind[,j])) > 10){
                stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders. \n
                     Please recode you confounder to satisfy this.")
            }
        }
    }

    geneffect = "additive"

    # Model variations

    # 1: Hardy-Weinberg equilibrium
    # 1.1: no interactions
    # 1.1.1: no xe,
    # 1.1.1.1: maf given
    # 1.1.1.2: maf not given
    # 1.1.2: xe given
    # 1.1.2.1: maf given
    # 1.1.2.2: maf not given
    # 1.2: gene-environment interactions
    # 1.2.1: maf given
    # 1.2.2: maf not given

    # 2: Hardy-Weinberg eqilibrium not assumed
    # 2.1: no interactions
    # 2.1.1: no xe,
    # 2.1.2: xe given
    # 2.1.2.1: no confounders
    # 2.1.2.2: confounders
    # 2.2: gene-environment interactions
    # 2.2.1: no confounders
    # 2.2.2: confounders

    if(hwe){
        #######################################################################
        # 1: Hardy-Weinberg
        #######################################################################
        if(interact == FALSE){
            ###################################################################
            # 1.1: No interactions
            ###################################################################
            if(!isxe){
                ###############################################################
                # 1.1.1: No xe
                ###############################################################
                if(!is.na(maf)){
                    ###########################################################
                    # 1.1.1.1: MAF given
                    ###########################################################
                    message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequency: ", toString(maf),sep = ""))
                    data = cbind(y,xg)
                    ng = dim(as.matrix(xg))[2]
                    model = epscomploglikmax(data, ng, hwe = TRUE, maf = maf,
                                             geneffect = "additive",
                                             hessian = TRUE)
                    hessian = model[[1]]
                    info = -1*ginv(hessian)
                    params = model[[2]]
                    sigma = params[length(params)]
                    coef = c()
                    ci = data.frame(matrix(NA,nrow = (length(params)-1),
                                           ncol = 2))
                    for(i in 1:(length(params)-1)){
                        coef[i] = params[i]
                        ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                        ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                    }
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)", newnames)
                    rownames(ci) = c("(intercept)",newnames)
                    result = list(coef,ci,sigma)
                    names(result) = c("coefficients","ci","sigma")
                    return(result)
                }else{
                    ###########################################################
                    # 1.1.1.2: MAF not given
                    ###########################################################
                    message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
                    data = cbind(y,xg)
                    ng = dim(as.matrix(xg))[2]
                    model = epscomploglikmax(data, ng, hwe = TRUE,
                                             geneffect = "additive",
                                             hessian = TRUE)
                    params = model[[2]]
                    nparam = 1 + ng
                    sigma = params[(nparam + 1)]
                    maf = params[(nparam + 2):length(params)]
                    hessian = model[[1]]
                    info = -1*ginv(hessian)
                    coef = c()
                    ci = data.frame(matrix(NA,nrow = (length(params)-1-ng),
                                           ncol = 2))
                    for(i in 1:nparam){
                        coef[i] = params[i]
                        ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                        ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                    }
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)",newnames)
                    rownames(ci) = c("(intercept)",newnames)
                    result = list(coef,ci,sigma,maf)
                    names(result) = c("coefficients","ci","sigma","maf")
                    return(result)
                }
            }else{
                ###############################################################
                # 1.1.2: xe
                ###############################################################
                data = cbind(y,xe,xg)
                xe = as.matrix(xe)
                xg = as.matrix(xg)
                ng = dim(as.matrix(xg))[2]
                ne = dim(as.matrix(xe))[2]

                if(!is.na(as.matrix(maf)[1,1])){
                    ###########################################################
                    # 1.1.2.1: MAF given
                    ###########################################################
                    message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequencies: ", toString(maf),sep = ""))
                    model = epscomploglikmax(data,ng, hwe = TRUE, maf = maf,
                                             geneffect = "additive",
                                             hessian = TRUE)
                    params = model[[2]]
                    nparam = 1 + ne + ng
                    sigma = params[(nparam + 1)]
                    hessian = model[[1]]
                    info = -1*ginv(hessian)
                    coef = c()
                    ci = data.frame(matrix(NA,nrow = (1 + ne + ng), ncol = 2))
                    for(i in 1:(1 + ne + ng)){
                        coef[i] = params[i]
                        ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                        ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                    }
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)",newnames)
                    rownames(ci) = c("(intercept)",newnames)

                    result = list(coef,ci,sigma)
                    names(result) = c("coefficients","ci","sigma")
                    return(result)
                }else{
                    ###########################################################
                    # 1.1.2.2: MAF not given
                    ###########################################################
                    message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
                    model = epscomploglikmax(data, ng, hwe = TRUE,
                                             geneffect = "additive",
                                             hessian = TRUE)
                    params = model[[2]]
                    nparam = 1 + ne + ng
                    sigma = params[(nparam + 1)]
                    maf = params[(nparam + 2):length(params)]
                    hessian = model[[1]]
                    info = -1*ginv(hessian)
                    coef = c()
                    ci = data.frame(matrix(NA,nrow = (length(params)-1-ng),
                                           ncol = 2))
                    for(i in 1:(1 + ne + ng)){
                        coef[i] = params[i]
                        ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                        ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                    }
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)",newnames)
                    rownames(ci) = c("(intercept)",newnames)
                    result = list(coef,ci,sigma,maf)
                    names(result) = c("coefficients","ci","sigma","maf")
                    return(result)
                }
            }
        }else{
            ###################################################################
            # 1.2: Interactions
            ###################################################################
            data = cbind(y,xe,xg)
            xe = as.matrix(xe)
            xg = as.matrix(xg)
            ng = dim(xg)[2]
            ne = dim(xe)[2]
            intnames = c()
            for(i in 1:length(interactind)){
                intnames[i] = paste(colnames(xe)[interactind[[i]][2]],"*",
                                    colnames(xg)[interactind[[i]][1]],sep="")
            }
            if(!is.na(as.matrix(maf)[1,1])){
                ###########################################################
                # 1.2.1: MAF given
                ###########################################################
                message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequency: ", toString(maf),sep = ""))
                model = epscomploglikmaxint(data, ng, interactind,
                                            hwe = TRUE, maf = maf,
                                            geneffect = "additive",
                                            hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ne + ng + length(interactind)
                sigma = params[(nparam + 1)]
                hessian = model[[1]]
                info = -1*ginv(hessian)
                coef = c()
                ci = data.frame(matrix(NA,nrow = (nparam), ncol = 2))
                for(i in 1:(nparam)){
                    coef[i] = params[i]
                    ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                    ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                }
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",newnames)
                rownames(ci) = c("(intercept)",newnames)
                result = list(coef,ci,sigma)
                names(result) = c("coefficients","ci","sigma")
                return(result)
            }else{
                ###########################################################
                # 1.2.1: MAF not given
                ###########################################################
                message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
                model = epscomploglikmaxint(data, ng, interactind,
                                            hwe = TRUE,
                                            geneffect = "additive",
                                            hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ne + ng + length(interactind)
                sigma = params[(nparam + 1)]
                maf = params[(nparam + 2):length(params)]
                hessian = model[[1]]
                info = -1*ginv(hessian)
                coef = c()
                ci = data.frame(matrix(NA,nrow = (nparam), ncol = 2))
                for(i in 1:(nparam)){
                    coef[i] = params[i]
                    ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                    ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                }
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",newnames)
                rownames(ci) = c("(intercept)",newnames)
                result = list(coef,ci,sigma,maf)
                names(result) = c("coefficients","ci","sigma","maf")
                return(result)
            }
        }
    }else{
        #####################################################################
        # 2. Hardy-Weinberg not assumed
        #####################################################################
        if(interact == FALSE){
            #################################################################
            # 2.1.1: No interactions
            #################################################################
            if(!isxe){
                ###############################################################
                # 2.1.1: No xe
                ###############################################################
                data = cbind(y,xg)
                ng = dim(as.matrix(xg))[2]
                model = epscomploglikmax(data, ng, gfreq = gfreq,
                                         geneffect = geneffect,
                                         hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ng
                sigma = params[(nparam + 1)]
                probs = params[(nparam + 2):length(params)]
                hessian = model[[1]]
                info = -1*ginv(hessian)
                coef = c()
                ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
                for(i in 1:nparam){
                    coef[i] = params[i]
                    ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
                    ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
                }
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",newnames)
                rownames(ci) = c("(intercept)",newnames)
                result = list(coef,ci,sigma,probs)
                names(result) = c("coefficients","ci","sigma","P(Xg)")
                return(result)
            }else{
                ###############################################################
                # 2.1.1: xe
                ###############################################################
                data = cbind(y,xe,xg)
                xe = as.matrix(xe)
                xg = as.matrix(xg)
                ng = dim(as.matrix(xg))[2]
                ne = dim(as.matrix(xe))[2]
                if(!confounder){
                    ###########################################################
                    # 2.1.1: No confounders
                    ###########################################################
                    model = epscomploglikmax(data, ng,
                                             gfreq = gfreq,
                                             geneffect = "additive",
                                             hessian = TRUE)

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
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)",
                                    newnames)
                    rownames(ci) = c("(intercept)",
                                     newnames)
                    result = list(coef,ci,sigma)
                    names(result) = c("coefficients","ci","sigma")
                    return(result)
                }else{
                    #########################################################
                    # 2.1.1: Confounders
                    #########################################################
                    model = epscomploglikmaxcond(data, ng,cind = cind,
                                                 hessian = TRUE)
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
                    colnames(ci) = c("lower 95% ci", "upper 95% ci")
                    names(coef) = c("(intercept)",
                                    newnames)
                    rownames(ci) = c("(intercept)",
                                     newnames)
                    result = list(coef,ci,sigma)
                    names(result) = c("coefficients","ci","sigma")
                    return(result)
                }
            }
        }else{
            #################################################################
            # 2.1.1: Interactions
            #################################################################
            data = cbind(y,xe,xg)
            xe = as.matrix(xe)
            xg = as.matrix(xg)
            ng = dim(as.matrix(xg))[2]
            ne = dim(as.matrix(xe))[2]
            intnames = c()
            for(i in 1:length(interactind)){
                intnames[i] = paste(colnames(xe)[interactind[[i]][2]],"*",
                                    colnames(xg)[interactind[[i]][1]],sep="")
            }
            if(!confounder){
                ###############################################################
                # 2.1.1: No confounders
                ###############################################################
                model = epscomploglikmaxint(data, ng, interactind = interactind,
                                            geneffect = geneffect, hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ne + ng + length(interactind)
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
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",
                                newnames)
                rownames(ci) = c("(intercept)",
                                 newnames)
                result = list(coef,ci,sigma)
                names(result) = c("coefficients","ci","sigma")
                return(result)
            }else{
                ###############################################################
                # 2.1.1: Confounders
                ###############################################################
                model = epscomploglikmaxcondint(data, ng,
                                                cind = cind,
                                                interactind = interactind,
                                                hessian = TRUE)
                params = model[[2]]
                nparam = 1 + ne + ng + length(interactind)
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
                colnames(ci) = c("lower 95% ci", "upper 95% ci")
                names(coef) = c("(intercept)",
                                newnames)
                rownames(ci) = c("(intercept)",
                                 newnames)
                result = list(coef,ci,sigma)
                names(result) = c("coefficients","ci","sigma")
                return(result)
            }
        }
    }
}
