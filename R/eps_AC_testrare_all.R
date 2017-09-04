
epsAC.rv.test.naive = function(epsdata0,covariates0,xg,isxe,confounder,cind){
    y = epsdata0[,1]
    n = length(y)
    xg = as.matrix(xg)

    if(isxe & !confounder){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # No confounding assumed
        ###############################################################
        xe = as.matrix(covariates0)
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = cbind(rep(1,n),xe)
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)

        gm = xg
        gm[!extreme,] = eg

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
        var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

        Sigma = var1 - var2

        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        t = t(s)%*%ginv(Sigma)%*%s
        pvalue = pchisq(t,ng,lower.tail=FALSE)

        result = list(t,pvalue)
        names(result) = c("statistic","p.value")
        return(result)

    }else if(isxe & confounder){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # Confounding assumed
        ###############################################################
        xe = as.matrix(covariates0)
        fit = lm(y~xe) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = cbind(rep(1,n),xe)
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

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

        varg = list()
        egmat =  matrix(0,ncol=ng,nrow=n)

        for(u in 1:nu){
            uind = uindex_cc[[u]]
            allind = uindex_all[[u]]
            for(c in 1:ng){
                egmat[allind,c] = mean(g[uind,c])
            }
            varg[[u]] = var(g[uind,])
        }

        gm = xg
        gm[!extreme,] = egmat[!extreme,]

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2

        var2 = 0*var1
        for(u in 1:nu){
            uind_ic = uindex_ic[[u]]
            uind_cc = uindex_cc[[u]]
            uind_all = uindex_all[[u]]

            n_uic = length(z[!extreme][uind_ic])
            n_ucc = length(z[extreme][uind_cc])

            var2 = var2 + varg[[u]]*(sum(z[!extreme][uind_ic]^2) - sigma2*(n_uic) - (1/n_ucc)*(sum(z[!extreme][uind_ic]))^2)/(sigma2^2)
        }

        Sigma = var1 - var2
        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        t = t(s)%*%ginv(Sigma)%*%s
        pvalue = pchisq(t,ng,lower.tail=FALSE)

        result = list(t,pvalue)
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

        Xe = rep(1,n)
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)

        gm = xg
        gm[!extreme,] = eg

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
        var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

        Sigma = var1 - var2

        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        t = t(s)%*%ginv(Sigma)%*%s
        pvalue = pchisq(t,ng,lower.tail=FALSE)

        result = list(t,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }
}

epsAC.rv.test.varcomp = function(epsdata0,covariates0,xg,isxe,confounder,cind,weights){

    Bmat = diag(weights)

    y = epsdata0[,1]
    n = length(y)
    xg = as.matrix(xg)
    xe = as.matrix(covariates0)

    if(isxe & !confounder){
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

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)

        gm = xg
        gm[!extreme,] = eg

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
        var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

        Sigma = var1 - var2

        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        #eg = colMeans(g)
        #varg = var(g)
        #egmat = matrix(eg,ncol=ng,nrow=n,byrow = TRUE)

        #sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        #sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

        #Sigma = sig1 + sig2
        #s = (z[extreme]%*%g +sum(z[!extreme])*eg)/sigma2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)

        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        accuracy = 1e6
        qdist = rep(0,accuracy)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(accuracy,1)
        }
        #critval = quantile(qdist,probs = 0.95)
        teststat = c(t(s)%*%Bmat%*%s)
        pvalue = sum(qdist > teststat)/accuracy


        result = list(teststat,pvalue,d)
        names(result) = c("statistic","p.value","eigen")
        return(result)
    }else if(isxe & confounder){
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

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

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

        varg = list()
        egmat =  matrix(0,ncol=ng,nrow=n)

        for(u in 1:nu){
            uind = uindex_cc[[u]]
            allind = uindex_all[[u]]
            for(c in 1:ng){
                egmat[allind,c] = mean(g[uind,c])
            }
            varg[[u]] = var(g[uind,])
        }

        gm = xg
        gm[!extreme,] = egmat[!extreme,]

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2

        var2 = 0*var1
        for(u in 1:nu){
            uind_ic = uindex_ic[[u]]
            uind_cc = uindex_cc[[u]]
            uind_all = uindex_all[[u]]

            n_uic = length(z[!extreme][uind_ic])
            n_ucc = length(z[extreme][uind_cc])

            var2 = var2 + varg[[u]]*(sum(z[!extreme][uind_ic]^2) - sigma2*(n_uic) - (1/n_ucc)*(sum(z[!extreme][uind_ic]))^2)/(sigma2^2)
        }

        Sigma = var1 - var2
        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        # extreme = !is.na(xg[,1])
        # g = as.matrix(xg[extreme,])
        # n_cc = sum(extreme)
        # ng = dim(g)[2]
        #
        # # confounding
        # ux = as.matrix(unique(xe[extreme,cind]))
        # nu = dim(ux)[1]
        # if(nu != dim(as.matrix(unique(xe[,cind])))[1]){
        #     warning("All unique levels of confounder not found in extreme sample")
        # }
        # uindex = list()
        # for(u in 1:nu){
        #     uindex[[u]] = (xe[extreme,cind] == ux[u,])
        # }
        # eg = list()
        # varg =  list()
        # egmat =  matrix(0,ncol=ng,nrow=n)
        # for(u in 1:nu){
        #     uind = uindex[[u]]
        #     allind = (xe[,cind] == ux[u,])
        #
        #     egmat[allind,] = colMeans(g[uind,])
        #
        #     eg[[u]] = colMeans(g[uind,])
        #     varg[[u]] = var(g[uind,])
        # }
        # sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        # sig2 = 0
        # s = z[extreme]%*%g
        # for(u in 1:nu){
        #     uind = uindex[[u]]
        #     n_uc = length(z[extreme][uind])
        #     sig2 = sig2 + varg[[u]]*(sum(z[extreme][uind]^2) - (1/n_uc)*(sum(z[extreme][uind]))^2)/(sigma2^2)
        #
        #     tempind = (xe[!extreme,cind] == ux[u,])
        #     s = s + sum(z[!extreme][tempind])*eg[[u]]
        # }
        # s = s/sigma2
        # Sigma = sig1 + sig2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)
        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        accuracy = 1e6
        qdist = rep(0,accuracy)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(accuracy,1)
        }
        #critval = quantile(qdist,probs = 0.95)
        teststat = c(t(s)%*%Bmat%*%s)
        pvalue = sum(qdist > teststat)/accuracy

        result = list(teststat,pvalue,d)
        names(result) = c("statistic","p.value","eigen")
        return(result)
    }else{
        ##############################################################
        # no environmental covariates (xe) present in the null model
        ##############################################################
        fit = lm(y~1) # Fit under H0
        z = fit$residuals
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = rep(1,n)
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)

        gm = xg
        gm[!extreme,] = eg

        var1 = (t(gm)%*%(Ie-He)%*%gm)/sigma2
        var2 = varg*(sum(z[!extreme]^2) - sigma2*(n-n_cc) - (1/n_cc)*(sum(z[!extreme]))^2)/(sigma2^2)

        Sigma = var1 - var2

        s = (t(gm)%*%(Ie-He)%*%y)/sigma2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)

        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        accuracy = 1e6
        qdist = rep(0,accuracy)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(accuracy,1)
        }
        #critval = quantile(qdist,probs = 0.95)
        teststat = c(t(s)%*%Bmat%*%s)
        pvalue = sum(qdist > teststat)/accuracy

        result = list(teststat,pvalue,d)
        names(result) = c("statistic","p.value","eigen")
        return(result)
    }
}
