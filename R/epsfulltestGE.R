#' @title Test for gene-environment associations under the EPS-full design
#' @description
#' \code{epsfull.testGE} performs a likelihood ratio test for gene-environment
#' interaction variables under the EPS-full design
#' @param nullmodel an object of class \code{\link[stats]{formula}}, that
#' describes the linear regression model under the null hypothesis
#' @param GE a list of interactions, with the colon-symbol used to denote
#' interaction, e.g. \code{c("xe1:xg1","xe1:xg2")}
#' @param confounder \code{TRUE} if distribution of SNPs should be
#' assumed dependent on other (non-genetic) covariates,
#' default set to \code{FALSE}
#' @param cx optional vector of names of confounding (non-genetic) covariates
#' @param hwe \code{TRUE} if Hardy-Weinberg equilibrium is assumed, default
#' set to \code{FALSE}
#' @param maf optional value for the minor allele frequencies under HWE
#' @param lrt use likelihood ratio test instead of score
#' @return \code{epsfull.testGE} returns
#' \item{statistic}{the value of the score test statistic}
#' \item{parameter}{the degrees of freedom of the statistic}
#' \item{p.value}{the P-value for the test}
#' @details
#' The \code{nullmodel} \code{\link[stats]{formula}} object is of the type
#' y~xe+xg, which describes a regression model, y=a+be*xe+bg*xg+e
#' assuming a normal distribution for the residuals (e). The covariate
#' xe is a non-genetic/environmental covariate (optional).
#' The covariate xg is a SNP (single-nucleotide polymorphism).
#' The variables are taken from the environment that the
#' function is called from.
#' Both xe and xg can be matrices.
#'
#' The test considers a regression model y=a+be*xe+bg*xg+b*xe*xg+e,
#' where b=0 under the null hypothesis. The output
#' of the function gives the test statistic and p-value for the test of
#' H0: b=0. The specific gene-environment interactions that should be
#' tested is specified in \code{GE}.
#'
#' The EPS-full design is such that the SNP genotype is only observed
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
#' xg1[!extreme] = NA
#' xg2[!extreme] = NA
#' xg = as.matrix(cbind(xg1,xg2))
#' xe = as.matrix(cbind(xe1,xe2))
#' epsfull.testGE(y~xe1+xe2+xg1+xg2,GE = c("xe1:xg1"))$p.value

epsAC.testGE = function(nullmodel, GE,confounder = FALSE, cx, hwe = FALSE, maf,lrt=FALSE){

    if(!lrt){
        #######################################################################
        # SCORE TEST
        #######################################################################

        if(class(nullmodel)!="formula"){
            stop("First argument must be of class formula")}

        options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
        options(na.action="na.omit")

        y = epsdata0[,1]
        n = length(y)

        covnames = colnames(covariates0)

        # which main effects are SNPs and which are environment?
        snpid = c()
        xid = c()
        ncov = length(covnames)
        for(c in 1:ncov){
            if(sum(is.na(covariates0[,c])) > 0){
                snpid[length(snpid)+1] = c
            }else if(sum(is.na(covariates0[,c])) == 0){
                xid[length(xid)+1] = c
            }
        }

        # which main effects should interact?
        interactind = list()
        nge = length(GE)
        xge = matrix(NA,ncol = nge,nrow = n)
        for(i in 1:nge){
            t = match(strsplit(GE[i],":")[[1]], covnames)
            if(sum(is.na(t)>0)){
                stop(paste("Cannot include interaction term ",toString(GE[i]),
                           ", because corresponding main effects were not found in the null model",
                           sep = ""))
            }
            xge[,i] = covariates0[,t[1]]*covariates0[,t[2]]
            if(t[1] %in% snpid){
                # always genetic variant first
                interactind[[i]] = c(match(t[1],snpid),match(t[2],xid))
            }else if(t[1] %in% xid){
                interactind[[i]] = c(match(t[2],snpid),match(t[1],xid))
            }
        }

        xg = as.matrix(covariates0[,snpid])
        xe = as.matrix(covariates0[,xid])
        colnames(xe) = colnames(covariates0)[xid]
        colnames(xg) = colnames(covariates0)[snpid]
        data = cbind(y,xe,xg)
        ng = dim(xg)[2]
        ne = dim(xe)[2]

        if(ng == 0){stop("No genetic covariates (with NAs) found")}


        if(confounder){
            if(missing(cx)){cx = colnames(xe)
            }else if(is.na(match(cx,colnames(xe)))){
                stop("The name of the confounder given is not the name of a covariate.")
            }
            message(paste("Confounders: ", toString(cx),sep=""))
            cind = match(cx,colnames(xe))
            xecind = as.matrix(xe[,cind])
            for(j in 1:length(cx)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                         Please recode you confounder to satisfy this.")
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

        genotypes = unique(g_cc)
        nug = dim(genotypes)[1]

        if(confounder){
            ux = as.matrix(unique(x_cc[,cind]))
            nu = dim(ux)[1]
            if(nu != dim(as.matrix(unique(xe[,cind])))[1]){
                warning("All unique levels of confounder not found in extreme sample")
            }
            uindex = list()
            for(u in 1:nu){
                uindex[[u]] = colSums((data.frame(t(x_ic[,cind])) == t(ux[u,])))
            }

            uindex_cc = list()
            for(u in 1:nu){
                uindex_cc[[u]] = colSums((data.frame(t(x_cc[,cind])) == t(ux[u,])))
            }
        }

            ###################################################################
            # Test GEIs one by one
            ###################################################################
            statistic = matrix(NA,ncol = 1, nrow = nge)
            parameter = matrix(NA,ncol = 1, nrow = nge)
            pvalue = matrix(NA,ncol = 1, nrow = nge)
            rownames(statistic) = GE
            rownames(parameter) = GE
            rownames(pvalue) = GE
            colnames(statistic) = "t"
            colnames(parameter) = "d.o.f"
            colnames(pvalue) = "p.value"
            if(hwe){
                ###############################################################
                # Hardy-Weinberg assumed, independent SNPs
                ###############################################################
            }else{
                ###############################################################
                # Hardy-Weinberg not assumed, SPNs can be dependent
                ###############################################################
                if(!confounder){
                    ###############################################################
                    # No confounding assumed
                    ###############################################################
                    fit0 = epsfullloglikmax(data,ng,hessian = TRUE)
                    coefs = fit0[[2]][1:length(fit0[[2]])-1]
                    alpha = coefs[1]
                    betaE = coefs[2:(2+ne-1)]
                    betaG = coefs[(2+ne):length(coefs)]
                    sigma = fit0[[2]][length(fit0[[2]])]
                    sigma2 = sigma*sigma

                    genoprobs = fit0[[3]]
                    genotypes = fit0[[4]]

                    for(i in 1:nge){

                        gint = interactind[[i]][1]
                        xint = interactind[[i]][2]
                        eg = g_cc[,gint]*x_cc[,xint]

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
                        I_14 = (1/sigma^2)*sum(g_cc[,gint]*x_cc[,xint]) +
                            (1/sigma^2)*t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs))%*%x_ic[,xint] -
                            (1/sigma^4)*t(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                              hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs)*
                                              hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%x_ic[,xint]

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
                        eg = g_cc[,gint]*x_cc[,xint]
                        I_24 = (1/sigma^2)*t(x_cc)%*%eg +
                            (1/sigma^2)*t(x_ic)%*%(x_ic[,xint]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs)) -
                            (1/sigma^4)*t(x_ic)%*%(x_ic[,xint]*(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
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
                        I_34 = (1/sigma^2)*t(eg)%*%g_cc + (1/sigma^2)*t(x_ic[,xint])%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs) -
                            (1/sigma^4)*(t(x_ic[,xint])%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs) -
                                             t(x_ic[,xint]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))%*%hfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs))
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
                            (1/sigma^2)*t(x_ic[,xint]^2)%*%hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,2,gint,genotypes,genoprobs) -
                            (1/sigma^4)*t(x_ic[,xint]^2)%*%(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,2,gint,genotypes,genoprobs) -
                                                                hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs)^2)

                        # betaEG-sigma
                        I_45 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%eg +
                            (2/sigma^3)*t(x_ic[,xint])%*%hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs) -
                            (1/sigma^5)*t(x_ic[,xint])%*%(hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,1,gint,genotypes,genoprobs) -
                                                              hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs)*
                                                              hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs))

                        # betaEG-genoprobs
                        I_46 = c()
                        for(p in 1:(dim(genotypes)[1]-1)){
                            I_46[p] = -(1/sigma^2)*t(x_ic[,xint])%*%(hfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,p) -
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

                        (I22old - t(I12old)%*%ginv(I11old)%*%I12old)

                        s = (sum((y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)*eg) +
                                 sum(x_ic[,xint]*hfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs)))/sigma2
                        t = (s*s)/Sigma
                        pval = pchisq(t,1,lower.tail=FALSE)

                        statistic[i,] = t
                        parameter[i,] = 1
                        pvalue[i,] = pval
                    }
                    result = list(statistic,parameter,pvalue)
                    names(result) = c("statistic","parameter","p.value")
                    return(result)
                }else if(confounder){
                    ###############################################################
                    # Confounding assumed
                    ###############################################################
                    fit0 = epsfullloglikmaxcond(data,ng,cind,hessian = TRUE)
                    coefs = fit0[[2]][1:length(fit0[[2]])-1]
                    alpha = coefs[1]
                    betaE = coefs[2:(2+ne-1)]
                    betaG = coefs[(2+ne):length(coefs)]
                    sigma = fit0[[2]][length(fit0[[2]])]
                    sigma2 = sigma*sigma

                    genoprobs = fit0[[3]]
                    genotypes = as.matrix(fit0[[3]][[1]][,2:(2+ng-1)])

                    for(i in 1:nge){

                        gint = interactind[[i]][1]
                        xint = interactind[[i]][2]
                        eg = g_cc[,gint]*x_cc[,xint]

                        # alpha-alpha
                        I_11 = (n/sigma^2) - (1/sigma^4)*sum(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                                                                 chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)^2)

                        # alpha-betaE
                        I_12 = (1/sigma^2)*colSums(xe)-
                            (1/sigma^4)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex) -
                                              (chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)^2))%*%x_ic
                        # alpha-betaG
                        I_13 = (1/sigma^2)*colSums(g_cc) +
                            (1/sigma^2)*colSums(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,uindex)) -
                            (1/sigma^4)*(colSums(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)) -
                                             t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

                        # alpha-betaEG
                        I_14 = (1/sigma^2)*sum(g_cc[,gint]*x_cc[,xint]) +
                            (1/sigma^2)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex))%*%x_ic[,xint] -
                            (1/sigma^4)*t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
                                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)*
                                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))%*%x_ic[,xint]

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
                        eg = g_cc[,gint]*x_cc[,xint]
                        I_24 = (1/sigma^2)*t(x_cc)%*%eg +
                            (1/sigma^2)*t(x_ic)%*%(x_ic[,xint]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex)) -
                            (1/sigma^4)*t(x_ic)%*%(x_ic[,xint]*(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
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
                        I_33 = (1/sigma^2)*t(g_cc)%*%g_cc + (1/sigma^2)*matrix(colSums(chfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,uindex)),nrow=ng,ncol=ng) -
                            (1/sigma^4)*(matrix(colSums(chfun2(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)),nrow=ng,ncol=ng) -
                                             t(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

                        # betaG-betaEG
                        I_34 = (1/sigma^2)*t(eg)%*%g_cc + (1/sigma^2)*t(x_ic[,xint])%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,0,1,gint,genotypes,genoprobs,uindex) -
                            (1/sigma^4)*(t(x_ic[,xint])%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,2,1,gint,genotypes,genoprobs,uindex) -
                                             t(x_ic[,xint]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))
                        # betaG-sigma
                        I_35 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%g_cc +
                            (2/sigma^3)*colSums(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)) -
                            (1/sigma^5)*(colSums(chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,3,0,gint,genotypes,genoprobs,uindex)) -
                                             t(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex))%*%
                                             chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex))

                        # betaG-genoprobs
                        I_36 = matrix(ncol=(nu*(dim(genotypes)[1]-1)),nrow=ng)
                        mati = 1
                        for(u in 1:nu){
                            t_36 = matrix(ncol=(dim(genotypes)[1]-1),nrow=ng)
                            for(p in 1:(dim(genotypes)[1]-1)){
                                t_36[,p] = -(1/sigma^2)*colSums(chfun1dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,p,uindex,u)) +
                                    (1/sigma^2)*t(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,0,0,gint,genotypes,genoprobs,p,uindex,u))%*%chfun1(y_ic,x_ic,alpha,betaE,betaG,sigma,1,0,gint,genotypes,genoprobs,uindex)
                            }
                            I_36[,mati:(mati+(dim(genotypes)[1]-1)-1)] = t_36
                            mati = mati + (dim(genotypes)[1]-1)
                        }


                        # betaEG-betaEG
                        I_44 = (1/sigma^2)*t(eg)%*%eg +
                            (1/sigma^2)*t(x_ic[,xint]^2)%*%chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,0,2,gint,genotypes,genoprobs,uindex) -
                            (1/sigma^4)*t(x_ic[,xint]^2)%*%(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,2,gint,genotypes,genoprobs,uindex) -
                                                                chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex)^2)

                        # betaEG-sigma
                        I_45 = (2/sigma^3)*t(y_cc - alpha - x_cc%*%betaE - g_cc%*%betaG)%*%eg +
                            (2/sigma^3)*t(x_ic[,xint])%*%chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex) -
                            (1/sigma^5)*t(x_ic[,xint])%*%(chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,3,1,gint,genotypes,genoprobs,uindex) -
                                                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,2,0,gint,genotypes,genoprobs,uindex)*
                                                              chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex))

                        # betaEG-genoprobs
                        I_46 = c()
                        for(u in 1:nu){
                            for(p in 1:(dim(genotypes)[1]-1)){
                                I_46[(length(I_46)+1)] = -(1/sigma^2)*t(x_ic[,xint])%*%(chfun0dash(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,p,uindex,u) -
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
                            cg_cc = as.matrix(g_cc[uindex_cc[[u]],])
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
                                 sum(x_ic[,xint]*chfun0(y_ic,x_ic,alpha,betaE,betaG,sigma,1,1,gint,genotypes,genoprobs,uindex)))/sigma2
                        t = (s*s)/Sigma
                        pval = pchisq(t,1,lower.tail=FALSE)

                        statistic[i,] = t
                        parameter[i,] = 1
                        pvalue[i,] = pval
                    }
                    result = list(statistic,parameter,pvalue)
                    names(result) = c("statistic","parameter","p.value")
                    return(result)
                }
        }

    }else{
        #######################################################################
        # Likelihood Ratio TEST
        #######################################################################
        onebyone = TRUE

        if(class(nullmodel)!="formula"){
            stop("First argument must be of class formula")}

        options(na.action="na.pass")
        epsdata0 = model.frame(nullmodel)
        covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
        options(na.action="na.omit")

        y = epsdata0[,1]
        n = length(y)

        covnames = colnames(covariates0)

        # which main effects are SNPs and which are environment?
        snpid = c()
        xid = c()
        ncov = length(covnames)
        for(c in 1:ncov){
            if(sum(is.na(covariates0[,c])) > 0){
                snpid[length(snpid)+1] = c
            }else if(sum(is.na(covariates0[,c])) == 0){
                xid[length(xid)+1] = c
            }
        }

        # which main effects should interact?
        interactind = list()
        nge = length(GE)
        xge = matrix(NA,ncol = nge,nrow = n)
        for(i in 1:nge){
            t = match(strsplit(GE[i],":")[[1]], covnames)
            if(sum(is.na(t)>0)){
                stop(paste("Cannot include interaction term ",toString(GE[i]),
                           ", because corresponding main effects were not found in the null model",
                           sep = ""))
            }
            xge[,i] = covariates0[,t[1]]*covariates0[,t[2]]
            if(t[1] %in% snpid){
                # always genetic variant first
                interactind[[i]] = c(match(t[1],snpid),match(t[2],xid))
            }else if(t[1] %in% xid){
                interactind[[i]] = c(match(t[2],snpid),match(t[1],xid))
            }
        }



        #if(confounder){stop("Confounding currently not allowed for interactions")}

        #geneffect = "additive"
        #if(missing(gfreq)){gfreq = NA}

        xg = as.matrix(covariates0[,snpid])
        xe = as.matrix(covariates0[,xid])
        colnames(xe) = colnames(covariates0)[xid]
        colnames(xg) = colnames(covariates0)[snpid]
        data = cbind(y,xe,xg)
        ng = dim(xg)[2]

        if(confounder){
            if(missing(cx)){cx = colnames(xe)
            }else if(is.na(match(cx,colnames(xe)))){
                stop("The name of the confounder given is not the name of a covariate.")
            }
            message(paste("Confounders: ", toString(cx),sep=""))
            cind = match(cx,colnames(xe))
            xecind = as.matrix(xe[,cind])
            for(j in 1:length(cx)){
                if(length(unique(xecind[,j])) > 10){
                    stop("Only discrete confounders with less than or equal to 10 unique levels are accepted as confounders.
                         Please recode you confounder to satisfy this.")
                }
                }
            }

        if(onebyone){
            #######################################################################
            # Test one interaction at a time
            #######################################################################
            nint = length(interactind)

            statistic = matrix(NA,ncol = nint, nrow = 1)
            parameter = matrix(NA,ncol = nint, nrow = 1)
            pvalue = matrix(NA,ncol = nint, nrow = 1)
            colnames(statistic) = GE
            colnames(parameter) = GE
            colnames(pvalue) = GE

            if(!confounder){
                if(hwe){
                    #######################################################
                    # Hardy Weinberg
                    #######################################################
                    if(missing(maf)){maf = NA}
                    fit0 = epsfullloglikmax(data,ng, hwe = TRUE, maf = maf,
                                            ll = TRUE)[[1]]
                    for(l in 1:nint){
                        fit1 = epsfullloglikmaxint(data,ng = ng,
                                                   hwe = TRUE, maf = maf,
                                                   interactind = list(interactind[[l]]),
                                                   ll = TRUE)[[1]]
                        t = -2*(fit0 - fit1)
                        pval = pchisq(t,1,lower.tail=FALSE)

                        statistic[,l] = t
                        parameter[,l] = 1
                        pvalue[,l] = pval
                    }
                    result = list(statistic,parameter,pvalue)
                    names(result) = c("statistic","parameter","p.value")
                    return(result)
                }else{
                    #######################################################
                    # Hardy Weinberg not assumed
                    #######################################################
                    fit0 = epsfullloglikmax(data,ng,
                                            ll = TRUE)[[1]]
                    for(l in 1:nint){
                        fit1 = epsfullloglikmaxint(data,ng = ng,
                                                   interactind = list(interactind[[l]]),
                                                   ll = TRUE)[[1]]
                        t = -2*(fit0 - fit1)
                        pval = pchisq(t,1,lower.tail=FALSE)

                        statistic[,l] = t
                        parameter[,l] = 1
                        pvalue[,l] = pval
                    }
                    result = list(statistic,parameter,pvalue)
                    names(result) = c("statistic","parameter","p.value")
                    return(result)
                }
            }else{
                fit0 = epsfullloglikmaxcond(data, ng,cind = cind,
                                            ll = TRUE)[[1]]
                for(l in 1:nint){
                    fit1 = epsfullloglikmaxcondint(data, ng,
                                                   cind = cind,
                                                   interactind = interactind,
                                                   ll = TRUE)[[1]]
                    t = -2*(fit0 - fit1)
                    pval = pchisq(t,1,lower.tail=FALSE)

                    statistic[,l] = t
                    parameter[,l] = 1
                    pvalue[,l] = pval
                }
                result = list(statistic,parameter,pvalue)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
        }else{
            #######################################################################
            # Test all interactions simultaneously
            #######################################################################
            if(hwe){
                #######################################################
                # Hardy Weinberg
                #######################################################
                if(missing(maf)){maf = NA}
                fit0 = epsfullloglikmax(data,ng, hwe = TRUE, maf = maf,
                                        ll = TRUE)[[1]]

                fit1 = epsfullloglikmaxint(data,ng = ng,
                                           hwe = TRUE, maf = maf,
                                           interactind = interactind,
                                           ll = TRUE)[[1]]
                t = -2*(fit0 - fit1)
                pval = pchisq(t,ng,lower.tail=FALSE)

                result = list(t,ng,pval)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }else{
                #######################################################
                # Hardy Weinberg not assumed
                #######################################################
                fit0 = epsfullloglikmax(data,ng,
                                        ll = TRUE)[[1]]

                fit1 = epsfullloglikmaxint(data,ng = ng,
                                           interactind = interactind,
                                           ll = TRUE)[[1]]
                t = -2*(fit0 - fit1)
                pval = pchisq(t,ng,lower.tail=FALSE)

                result = list(t,ng,pval)
                names(result) = c("statistic","parameter","p.value")
                return(result)
            }
        }
    }

}

