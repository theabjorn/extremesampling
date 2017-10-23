#' @title Test for gene-environment associations under the EPS-AC design
#' @description
#' \code{epsAC.test.GE} performs a likelihood ratio test for gene-environment
#' interaction variables under the EPS-AC design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' (including the genetic covariate)
#' @param G The genetic covariate in the GE interaction
#' @param E The environmental covariate in the GE interaction
#' @param confounder (optional) vector of names of confounding (non-genetic) covariates
#' @return \code{epsAC.testGE} returns
#' \item{statistic}{thescore test statistic}
#' \item{p.value}{the P-value}
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe+xg, which describes a regression model, y=a+be*xe+bg*xg+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The covariate xg is one or more, genetic markers where missing values
#' are coded as NA. Missing-mechanism must be MCAR or MAR.
#'
#' Confounders are discrete covariates (xe) and the distribution of xg is
#' modelled for each level of unique value of xe.
#'
#' The test considers a regression model y=a+be*xe+bg*xg+b*xe*xg+e,
#' where b=0 under the null hypothesis. The output
#' of the function gives the test statistic and p-value for the test of
#' H0: b=0. The genetic covarite (G) and the environmental covariate (E) for
#' the interaction must be specified, G with missing values coded as NA and
#' E with no missing values. G must also be included in the null model.
#'
#' @import MASS stats
#' @export
#' @examples
#'
#'

epsAC.test.GE = function(nullmodel,G,E,confounder){

    if(class(nullmodel)!="formula"){
        stop("First argument must be of class formula")}

    options(na.action="na.pass")
    epsdata0 = model.frame(nullmodel)
    covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
    options(na.action="na.omit")

    y = epsdata0[,1]
    n = length(y)

    G = as.matrix(G)
    E = as.matrix(E)
    if((dim(G)[2]>1)|(dim(E)[2]>1)){
        stop("G and E must be vectors of one genetic covariate and one environmental covariate")
    }

    covnames = colnames(covariates0)

    conf = FALSE
    if(!missing(confounder)){conf = TRUE}


    # which main effects are xgs and which are environment?
    xgid = c()
    xid = c()
    ncov = length(covnames)
    for(c in 1:ncov){
        if(sum(is.na(covariates0[,c])) > 0){
            xgid[length(xgid)+1] = c
        }else if(sum(is.na(covariates0[,c])) == 0){
            xid[length(xid)+1] = c
        }
    }

    xg = as.matrix(covariates0[,xgid])
    xe = as.matrix(covariates0[,xid])
    colnames(xe) = colnames(covariates0)[xid]
    colnames(xg) = colnames(covariates0)[xgid]
    data = cbind(y,xe,xg)
    ng = dim(xg)[2]
    ne = dim(xe)[2]

    if(ng == 0){stop("Genetic covariates must be included in null model")}
    if(ng > 1){stop("Only G should be included in the null model, no other genetic covariates")}

    isxe = TRUE
    if(ne == 0){isxe = FALSE}

    if(conf){
        if(is.na(match(confounder,colnames(xe)))){
            stop("The name of the confounder given is not the name of a covariate.")
        }
        message(paste("Confounders: ", toString(confounder),sep=""))
        cind = match(confounder,colnames(xe))
        xecind = as.matrix(xe[,cind])
        for(j in 1:length(confounder)){
            if(length(unique(xecind[,j])) > 10){
                stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.")
            }
        }
    }

    # Assume missing is the same for all xg
    extreme = !is.na(xg[,1])
    y_cc = y[extreme]
    y_ic = y[!extreme]
    g_cc = as.matrix(xg[extreme,])
    x_cc = as.matrix(xe[extreme,])
    x_ic = as.matrix(xe[!extreme,])
    n_cc = length(y_cc)
    n_ic = length(y_ic)

    if(conf){
        ux = as.matrix(unique(x_ic[,cind]))
        nu = dim(ux)[1]
        if(nu != dim(as.matrix(unique(xe[,cind])))[1]){
            warning("All unique levels of confounder not found in extreme sample")
        }
        uindex = list()
        uindex_cc = list()
        for(u in 1:nu){
            uindex[[u]] = colSums((data.frame(t(x_ic[,cind])) == t(ux[u,])))
            uindex_cc[[u]] = colSums((data.frame(t(x_cc[,cind])) == t(ux[u,])))
        }
    }


    if(!conf){
        ###############################################################
        # No confounding assumed
        ###############################################################
        if(isxe){
            fit0 = epsAC.loglikmax(data,ng,hessian = TRUE)
            coefs = fit0[[2]][1:length(fit0[[2]])-1]
            alpha = coefs[1]
            betaE = coefs[2:(2+ne-1)]
            betaG = coefs[(2+ne):length(coefs)]
            sigma = fit0[[2]][length(fit0[[2]])]
            sigma2 = sigma*sigma

            genoprobs = fit0[[3]]
            genotypes = fit0[[4]]

            eg = G[extreme,1]*E[extreme,1]
            gint = 1

            # alpha-alpha
            I_11 = (n/sigma^2) - (1/sigma^4)*sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs) -
                                                     hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)^2)

            # alpha-betaE
            I_12 = (1/sigma^2)*colSums(xe)-
                (1/sigma^4)*t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs) -
                                  (hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)^2))%*%x_ic
            # alpha-betaG
            I_13 = (1/sigma^2)*colSums(g_cc) +
                (1/sigma^2)*colSums(hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs)) -
                (1/sigma^4)*(colSums(hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)) -
                                 t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # alpha-betaEG
            I_14 = (1/sigma^2)*sum(G[extreme,]*E[extreme,]) +
                (1/sigma^2)*t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs))%*%E[!extreme,] -
                (1/sigma^4)*t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                  hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)*
                                  hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%E[!extreme,]

            # alpha-sigma
            I_15 = (2/sigma^3)*sum(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG) +
                (2/sigma^3)*sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)) -
                (1/sigma^5)*(sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs)) -
                                 sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                         hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)))

            # alpha-genoprobs
            I_16 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_16[p] = -(1/sigma^2)*sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^2)*sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))
            }

            # betaE-betaE
            I_22 = (1/sigma^2)*t(xe)%*%xe-
                (1/sigma^4)*t(x_ic)%*%diag(c(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs) -
                                                 (hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)^2)))%*%x_ic
            # betaE-betaG
            I_23 = (1/sigma^2)*t(x_cc)%*%g_cc +
                (1/sigma^2)*t(x_ic)%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs) -
                (1/sigma^4)*(t(x_ic)%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs) -
                                 t(x_ic)%*%diag(c(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # betaE-betaEG
            eg = G[extreme,]*E[extreme,]
            I_24 = (1/sigma^2)*t(x_cc)%*%eg +
                (1/sigma^2)*t(x_ic)%*%(E[!extreme,]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs)) -
                (1/sigma^4)*t(x_ic)%*%(E[!extreme,]*(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                                         hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs)))

            # betaE-sigma
            I_25 = (2/sigma^3)*t(x_cc)%*%c(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG) +
                (2/sigma^3)*t(x_ic)%*%c(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)) -
                (1/sigma^5)*(t(x_ic)%*%c(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs) -
                                             c(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                                   hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))))

            # betaE-genoprobs
            I_26 = matrix(ncol=(dim(genotypes)[1]-1),nrow=dim(xe)[2])
            for(p in 1:(dim(genotypes)[1]-1)){
                I_26[,p] = -(1/sigma^2)*t(x_ic)%*%(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p) -
                                                       hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))
            }

            # betaG-betaG
            I_33 = (1/sigma^2)*t(g_cc)%*%g_cc + (1/sigma^2)*matrix(colSums(hfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs)),nrow=ng,ncol=ng) -
                (1/sigma^4)*(matrix(colSums(hfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)),nrow=ng,ncol=ng) -
                                 t(hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # betaG-betaEG
            I_34 = (1/sigma^2)*t(eg)%*%g_cc + (1/sigma^2)*t(E[!extreme,])%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs) -
                (1/sigma^4)*(t(E[!extreme,])%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                 t(E[!extreme,]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))
            # betaG-sigma
            I_35 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%g_cc +
                (2/sigma^3)*colSums(hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)) -
                (1/sigma^5)*(colSums(hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs)) -
                                 t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs))%*%
                                 hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # betaG-genoprobs
            I_36 = matrix(ncol=(dim(genotypes)[1]-1),nrow=ng)
            for(p in 1:(dim(genotypes)[1]-1)){
                I_36[,p] = -(1/sigma^2)*colSums(hfun1dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^2)*t(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)
            }

            # betaEG-betaEG
            I_44 = (1/sigma^2)*t(eg)%*%eg +
                (1/sigma^2)*t(E[!extreme,]^2)%*%hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,2,gint,genotypes,genoprobs) -
                (1/sigma^4)*t(E[!extreme,]^2)%*%(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,2,gint,genotypes,genoprobs) -
                                                     hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs)^2)

            # betaEG-sigma
            I_45 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%eg +
                (2/sigma^3)*t(E[!extreme,])%*%hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs) -
                (1/sigma^5)*t(E[!extreme,])%*%(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,1,gint,genotypes,genoprobs) -
                                                   hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                                   hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))

            # betaEG-genoprobs
            I_46 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_46[p] = -(1/sigma^2)*t(E[!extreme,])%*%(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,p) -
                                                              hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))
            }

            # sigma-sigma

            I_55 = -n/sigma2 + (3/sigma^4)*sum((y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2) +
                (3/sigma^4)*sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)) -
                (1/sigma^6)*sum(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,4,0,gint,genotypes,genoprobs) -
                                    hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)^2)

            # sigma-genoprobs
            I_56 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_56[p] = -(1/sigma^3)*sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^3)*sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs))
            }

            # genoprobs-genoprobs
            mp = dim(genotypes)[1]
            I_66 = matrix(sum(as.vector(rowSums(g_cc - matrix(genotypes[mp,],ncol = dim(g_cc)[2],nrow=dim(g_cc)[1])))==0)/(genoprobs[mp]^2),ncol=(dim(genotypes)[1]-1),nrow=(dim(genotypes)[1]-1))
            #I_66 = matrix(0,ncol=(dim(genotypes)[1]-1),nrow=(dim(genotypes)[1]-1))
            for(p in 1:(dim(genotypes)[1]-1)){
                for(p2 in 1:(dim(genotypes)[1]-1)){
                    if(p==p2){
                        I_66[p,p] = I_66[p,p] + sum(as.vector(rowSums(g_cc - matrix(genotypes[p,],ncol = dim(g_cc)[2],nrow=dim(g_cc)[1])))==0)/(genoprobs[p]^2) +
                            sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)^2)
                    }else{
                        I_66[p,p2] = I_66[p,p2] + sum(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p2))
                    }
                }
            }

            I11_1 = c(I_11,I_12,I_13,I_15,I_16)
            I11_2 = cbind(t(I_12),I_22,I_23,I_25,I_26)
            I11_3 = cbind(t(I_13),t(I_23),I_33,t(I_35),I_36)
            I11_5 = c(I_15,t(I_25),I_35,I_55,I_56)
            I11_6 = cbind(I_16,t(I_26),t(I_36),I_56,I_66)
            I11 = rbind(I11_1,I11_2,I11_3,I11_5,I11_6)

            I22 = I_44

            I12 = c(I_14,t(I_24),I_34,I_45,I_46)

            Sigma = (I22 - t(I12)%*%ginv(I11)%*%I12)

            s = (sum((y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)*eg) +
                     sum(E[!extreme,]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs)))/sigma2
            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            result = list(t,pval)
            names(result) = c("statistic","p.value")
            return(result)
        }else{
            fit0 = epsAC.loglikmax(data,ng,hessian = TRUE)
            coefs = fit0[[2]][1:length(fit0[[2]])-1]
            alpha = coefs[1]
            betaG = coefs[2]
            sigma = fit0[[2]][length(fit0[[2]])]
            sigma2 = sigma*sigma

            genoprobs = fit0[[3]]
            genotypes = fit0[[4]]

            eg = G[extreme,1]*E[extreme,1]
            gint = 1
            x_ic = as.matrix(E[!extreme,])

            # alpha-alpha
            I_11 = (n/sigma^2) - (1/sigma^4)*sum(hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs) -
                                                     hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)^2)

            # alpha-betaG
            I_13 = (1/sigma^2)*sum(g_cc) +
                (1/sigma^2)*sum(hfun10(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs)) -
                (1/sigma^4)*(sum(hfun10(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)) -
                                 t(hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))%*%hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # alpha-betaEG
            I_14 = (1/sigma^2)*sum(G[extreme,]*E[extreme,]) +
                (1/sigma^2)*t(hfun00(y_ic,alpha,betaG,sigma,0,1,gint,genotypes,genoprobs))%*%E[!extreme,] -
                (1/sigma^4)*t(hfun00(y_ic,alpha,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                  hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)*
                                  hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%E[!extreme,]

            # alpha-sigma
            I_15 = (2/sigma^3)*sum(y_cc - alpha  - g_cc%*%betaG) +
                (2/sigma^3)*sum(hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)) -
                (1/sigma^5)*(sum(hfun00(y_ic,alpha,betaG,sigma,3,0,gint,genotypes,genoprobs)) -
                                 sum(hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                         hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)))

            # alpha-genoprobs
            I_16 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_16[p] = -(1/sigma^2)*sum(hfun0dash0(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^2)*sum(hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun00(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))
            }

            # betaG-betaG
            I_33 = (1/sigma^2)*t(g_cc)%*%g_cc + (1/sigma^2)*matrix(sum(hfun20(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs)),nrow=ng,ncol=ng) -
                (1/sigma^4)*(matrix(sum(hfun20(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)),nrow=ng,ncol=ng) -
                                 t(hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))%*%hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # betaG-betaEG
            I_34 = (1/sigma^2)*t(eg)%*%g_cc + (1/sigma^2)*t(E[!extreme,])%*%hfun10(y_ic,alpha,betaG,sigma,0,1,gint,genotypes,genoprobs) -
                (1/sigma^4)*(t(E[!extreme,])%*%hfun10(y_ic,alpha,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                 t(E[!extreme,]*hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))
            # betaG-sigma
            I_35 = (2/sigma^3)*t(y_cc - alpha  - g_cc%*%betaG)%*%g_cc +
                (2/sigma^3)*sum(hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)) -
                (1/sigma^5)*(sum(hfun10(y_ic,alpha,betaG,sigma,3,0,gint,genotypes,genoprobs)) -
                                 t(hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs))%*%
                                 hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs))

            # betaG-genoprobs
            I_36 = matrix(ncol=(dim(genotypes)[1]-1),nrow=ng)
            for(p in 1:(dim(genotypes)[1]-1)){
                I_36[,p] = -(1/sigma^2)*sum(hfun1dash0(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^2)*t(hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p))%*%hfun10(y_ic,alpha,betaG,sigma,1,0,gint,genotypes,genoprobs)
            }

            # betaEG-betaEG
            I_44 = (1/sigma^2)*t(eg)%*%eg +
                (1/sigma^2)*t(E[!extreme,]^2)%*%hfun00(y_ic,alpha,betaG,sigma,0,2,gint,genotypes,genoprobs) -
                (1/sigma^4)*t(E[!extreme,]^2)%*%(hfun00(y_ic,alpha,betaG,sigma,2,2,gint,genotypes,genoprobs) -
                                                     hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs)^2)

            # betaEG-sigma
            I_45 = (2/sigma^3)*t(y_cc - alpha  - g_cc%*%betaG)%*%eg +
                (2/sigma^3)*t(E[!extreme,])%*%hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs) -
                (1/sigma^5)*t(E[!extreme,])%*%(hfun00(y_ic,alpha,betaG,sigma,3,1,gint,genotypes,genoprobs) -
                                                   hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                                   hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs))

            # betaEG-genoprobs
            I_46 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_46[p] = -(1/sigma^2)*t(E[!extreme,])%*%(hfun0dash0(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs,p) -
                                                              hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs))
            }

            # sigma-sigma

            I_55 = -n/sigma2 + (3/sigma^4)*sum((y_cc - alpha  - g_cc%*%betaG)^2) +
                (3/sigma^4)*sum(hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)) -
                (1/sigma^6)*sum(hfun00(y_ic,alpha,betaG,sigma,4,0,gint,genotypes,genoprobs) -
                                    hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs)^2)

            # sigma-genoprobs
            I_56 = c()
            for(p in 1:(dim(genotypes)[1]-1)){
                I_56[p] = -(1/sigma^3)*sum(hfun0dash0(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs,p)) +
                    (1/sigma^3)*sum(hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun00(y_ic,alpha,betaG,sigma,2,0,gint,genotypes,genoprobs))
            }

            # genoprobs-genoprobs
            mp = dim(genotypes)[1]
            I_66 = matrix(sum(as.vector(rowSums(g_cc - matrix(genotypes[mp,],ncol = dim(g_cc)[2],nrow=dim(g_cc)[1])))==0)/(genoprobs[mp]^2),ncol=(dim(genotypes)[1]-1),nrow=(dim(genotypes)[1]-1))
            #I_66 = matrix(0,ncol=(dim(genotypes)[1]-1),nrow=(dim(genotypes)[1]-1))
            for(p in 1:(dim(genotypes)[1]-1)){
                for(p2 in 1:(dim(genotypes)[1]-1)){
                    if(p==p2){
                        I_66[p,p] = I_66[p,p] + sum(as.vector(rowSums(g_cc - matrix(genotypes[p,],ncol = dim(g_cc)[2],nrow=dim(g_cc)[1])))==0)/(genoprobs[p]^2) +
                            sum(hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p)^2)
                    }else{
                        I_66[p,p2] = I_66[p,p2] + sum(hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p)*hfun0dash0(y_ic,alpha,betaG,sigma,0,0,gint,genotypes,genoprobs,p2))
                    }
                }
            }

            I11_1 = c(I_11,I_13,I_15,I_16)
            I11_3 = cbind(t(I_13),I_33,t(I_35),I_36)
            I11_5 = c(I_15,I_35,I_55,I_56)
            I11_6 = cbind(I_16,t(I_36),I_56,I_66)
            I11 = rbind(I11_1,I11_3,I11_5,I11_6)

            I22 = I_44

            I12 = c(I_14,I_34,I_45,I_46)

            Sigma = (I22 - t(I12)%*%ginv(I11)%*%I12)

            s = (sum((y_cc - alpha - g_cc%*%betaG)*eg) +
                     sum(E[!extreme,]*hfun00(y_ic,alpha,betaG,sigma,1,1,gint,genotypes,genoprobs)))/sigma2
            t = (s*s)/Sigma
            pval = pchisq(t,1,lower.tail=FALSE)

            result = list(t,pval)
            names(result) = c("statistic","p.value")
            return(result)
        }

    }else if(conf){
        ###############################################################
        # Confounding assumed
        ###############################################################
        fit0 = epsAC.loglikmaxcond(data,ng,cind,hessian = TRUE)
        coefs = fit0[[2]][1:length(fit0[[2]])-1]
        alpha = coefs[1]
        betaE = coefs[2:(2+ne-1)]
        betaG = coefs[(2+ne):length(coefs)]
        sigma = fit0[[2]][length(fit0[[2]])]
        sigma2 = sigma*sigma

        genoprobs = fit0[[3]]
        genotypes = as.matrix(fit0[[3]][[1]][,2:(2+ng-1)])

        gint = 1
        eg = G[extreme,]*E[extreme,]


        # alpha-alpha
        I_11 = (n/sigma^2) - (1/sigma^4)*sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                                                 chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)^2)

        # alpha-betaE
        I_12 = (1/sigma^2)*sum(xe)-
            (1/sigma^4)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                              (chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)^2))%*%x_ic
        # alpha-betaG
        I_13 = (1/sigma^2)*sum(g_cc) +
            (1/sigma^2)*sum(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^4)*(sum(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)) -
                             t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

        # alpha-betaEG
        I_14 = (1/sigma^2)*sum(G[extreme,]*E[extreme,]) +
            (1/sigma^2)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex))%*%E[!extreme,] -
            (1/sigma^4)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)*
                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))%*%E[!extreme,]

        # alpha-sigma
        I_15 = (2/sigma^3)*sum(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG) +
            (2/sigma^3)*sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^5)*(sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs,uindex)) -
                             sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)*
                                     chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)))

        # alpha-genoprobs
        I_16 = c()
        for(u in 1:nu){
            for(p in 1:(dim(genotypes)[1]-1)){
                I_16[(length(I_16)+1)] = -(1/sigma^2)*sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p,uindex,u)) +
                    (1/sigma^2)*sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))
            }
        }

        # betaE-betaE
        I_22 = (1/sigma^2)*t(xe)%*%xe-
            (1/sigma^4)*t(x_ic)%*%diag(c(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                                             (chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)^2)))%*%x_ic
        # betaE-betaG
        I_23 = (1/sigma^2)*t(x_cc)%*%g_cc +
            (1/sigma^2)*t(x_ic)%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,uindex) -
            (1/sigma^4)*(t(x_ic)%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                             t(x_ic)%*%diag(c(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

        # betaE-betaEG
        eg = G[extreme,]*E[extreme,]
        I_24 = (1/sigma^2)*t(x_cc)%*%eg +
            (1/sigma^2)*t(x_ic)%*%(E[!extreme,]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^4)*t(x_ic)%*%(E[!extreme,]*(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
                                                     chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex)))

        # betaE-sigma
        I_25 = (2/sigma^3)*t(x_cc)%*%c(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG) +
            (2/sigma^3)*t(x_ic)%*%c(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^5)*(t(x_ic)%*%c(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs,uindex) -
                                         c(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)*
                                               chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))))

        # betaE-genoprobs
        I_26 = matrix(ncol=(nu*(dim(genotypes)[1]-1)),nrow=dim(xe)[2])
        counti = 0
        for(u in 1:nu){
            for(p in 1:(dim(genotypes)[1]-1)){
                counti = counti + 1
                I_26[,counti] = -(1/sigma^2)*t(x_ic)%*%(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p,uindex,u) -
                                                            chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))
            }
        }


        # betaG-betaG
        I_33 = (1/sigma^2)*t(g_cc)%*%g_cc + (1/sigma^2)*matrix(sum(chfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,uindex)),nrow=ng,ncol=ng) -
            (1/sigma^4)*(matrix(sum(chfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)),nrow=ng,ncol=ng) -
                             t(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

        # betaG-betaEG
        I_34 = (1/sigma^2)*t(eg)%*%g_cc + (1/sigma^2)*t(E[!extreme,])%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex) -
            (1/sigma^4)*(t(E[!extreme,])%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
                             t(E[!extreme,]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))
        # betaG-sigma
        I_35 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%g_cc +
            (2/sigma^3)*sum(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^5)*(sum(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs,uindex)) -
                             t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex))%*%
                             chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

        # betaG-genoprobs
        I_36 = matrix(ncol=(nu*(dim(genotypes)[1]-1)),nrow=ng)
        mati = 1
        for(u in 1:nu){
            t_36 = matrix(ncol=(dim(genotypes)[1]-1),nrow=ng)
            for(p in 1:(dim(genotypes)[1]-1)){
                t_36[,p] = -(1/sigma^2)*sum(chfun1dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p,uindex,u)) +
                    (1/sigma^2)*t(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)
            }
            I_36[,mati:(mati+(dim(genotypes)[1]-1)-1)] = t_36
            mati = mati + (dim(genotypes)[1]-1)
        }


        # betaEG-betaEG
        I_44 = (1/sigma^2)*t(eg)%*%eg +
            (1/sigma^2)*t(E[!extreme,]^2)%*%chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,2,gint,genotypes,genoprobs,uindex) -
            (1/sigma^4)*t(E[!extreme,]^2)%*%(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,2,gint,genotypes,genoprobs,uindex) -
                                                 chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex)^2)

        # betaEG-sigma
        I_45 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%eg +
            (2/sigma^3)*t(E[!extreme,])%*%chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex) -
            (1/sigma^5)*t(E[!extreme,])%*%(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,1,gint,genotypes,genoprobs,uindex) -
                                               chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)*
                                               chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))

        # betaEG-genoprobs
        I_46 = c()
        for(u in 1:nu){
            for(p in 1:(dim(genotypes)[1]-1)){
                I_46[(length(I_46)+1)] = -(1/sigma^2)*t(E[!extreme,])%*%(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,p,uindex,u) -
                                                                             chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))
            }
        }

        # sigma-sigma

        I_55 = -n/sigma2 + (3/sigma^4)*sum((y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)^2) +
            (3/sigma^4)*sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)) -
            (1/sigma^6)*sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,4,0,gint,genotypes,genoprobs,uindex) -
                                chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)^2)

        # sigma-genoprobs
        I_56 = c()
        for(u in 1:nu){
            for(p in 1:(dim(genotypes)[1]-1)){
                I_56[(length(I_56)+1)] = -(1/sigma^3)*sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,p,uindex,u)) +
                    (1/sigma^3)*sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex))
            }
        }

        # genoprobs-genoprobs
        mp = length(genotypes)
        I_66 = matrix(0,ncol=nu*(dim(genotypes)[1]-1),nrow = nu*(dim(genotypes)[1]-1))
        mati = 1
        for(u in 1:nu){
            cg_cc = as.matrix(g_cc[(uindex_cc[[u]]==1),])
            cgenoprobs = genoprobs[[u]][,1]
            t_66 = matrix(sum(as.vector(rowSums(cg_cc - matrix(genotypes[mp,],ncol = dim(cg_cc)[2],nrow=dim(cg_cc)[1])))==0)/(cgenoprobs[mp]^2),ncol=(dim(genotypes)[1]-1),nrow=(dim(genotypes)[1]-1))

            for(p in 1:(dim(genotypes)[1]-1)){
                for(p2 in 1:(dim(genotypes)[1]-1)){
                    if(p==p2){
                        t_66[p,p] = t_66[p,p] + sum(as.vector(rowSums(cg_cc - matrix(genotypes[p,],ncol = dim(cg_cc)[2],nrow=dim(cg_cc)[1])))==0)/(cgenoprobs[p]^2) +
                            sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)^2)
                    }else{
                        t_66[p,p2] = t_66[p,p2] + sum(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u)*chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p2,uindex,u))
                    }
                }
            }
            I_66[mati:(mati+(dim(genotypes)[1]-1)-1),mati:(mati+(dim(genotypes)[1]-1)-1)] = t_66
            mati = mati + (dim(genotypes)[1]-1)
        }

        I11_1 = c(I_11,I_12,I_13,I_15,I_16)
        I11_2 = cbind(t(I_12),I_22,I_23,I_25,I_26)
        I11_3 = cbind(t(I_13),t(I_23),I_33,t(I_35),I_36)
        I11_5 = c(I_15,t(I_25),I_35,I_55,I_56)
        I11_6 = cbind(I_16,t(I_26),t(I_36),I_56,I_66)
        I11 = rbind(I11_1,I11_2,I11_3,I11_5,I11_6)

        I22 = I_44

        I12 = c(I_14,t(I_24),I_34,I_45,I_46)

        Sigma = (I22 - t(I12)%*%ginv(I11)%*%I12)

        s = (sum((y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)*eg) +
                 sum(E[!extreme,]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex)))/sigma2
        t = (s*s)/Sigma
        pval = pchisq(t,1,lower.tail=FALSE)

        result = list(t,pval)
        names(result) = c("statistic","p.value")
        return(result)
    }
}
