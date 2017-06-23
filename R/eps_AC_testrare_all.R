
epsAC.rv.test.naive = function(epsdata0,covariates0,RV,isxe,confounder,cind){
    y = epsdata0[,1]
    n = length(y)
    xg = as.matrix(RV)

    if(isxe & !confounder){
        ###############################################################
        # Environmental covariates (xe) present in the null model
        # No confounding assumed
        ###############################################################
        xe = as.matrix(covariates0)
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

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)
        egmat = matrix(eg,ncol=ng,nrow=n,byrow = TRUE)

        sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

        Sigma = sig1 + sig2
        s = (z[extreme]%*%g +sum(z[!extreme])*eg)/sigma2

        t = s%*%ginv(Sigma)%*%t(s)
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
        coefs = fit$coef
        alpha = coefs[1]
        beta = coefs[2:length(coefs)]
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
        uindex = list()
        for(u in 1:nu){
            uindex[[u]] = (xe[extreme,cind] == ux[u,])
        }

        eg = list()
        varg =  list()
        egmat =  matrix(0,ncol=ng,nrow=n)

        for(u in 1:nu){
            uind = uindex[[u]]
            allind = (xe[,cind] == ux[u,])

            egmat[allind,] = colMeans(g[uind,])

            eg[[u]] = colMeans(g[uind,])
            varg[[u]] = var(g[uind,])
        }

        sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        sig2 = 0

        s = z[extreme]%*%g

        for(u in 1:nu){
            uind = uindex[[u]]
            n_uc = length(z[extreme][uind])
            sig2 = sig2 + varg[[u]]*(sum(z[extreme][uind]^2) - (1/n_uc)*(sum(z[extreme][uind]))^2)/(sigma2^2)

            tempind = (xe[!extreme,cind] == ux[u,])
            s = s + sum(z[!extreme][tempind])*eg[[u]]
        }

        s = s/sigma2
        Sigma = sig1 + sig2

        t = s%*%ginv(Sigma)%*%t(s)
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
        coefs = fit$coef
        alpha = coefs[1]
        sigma2 = sum(z^2)/n
        sigma = sqrt(sigma2)

        Xe = cbind(rep(1,n))
        He = Xe%*%ginv(t(Xe)%*%Xe)%*%t(Xe)
        Ie = diag(1,nrow = n, ncol = n)

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)
        egmat = matrix(eg,ncol=ng,nrow=n,byrow = TRUE)

        sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

        Sigma = sig1 + sig2
        s = (z[extreme]%*%g +sum(z[!extreme])*eg)/sigma2

        t = s%*%ginv(Sigma)%*%t(s)
        pvalue = pchisq(t,ng,lower.tail=FALSE)

        result = list(t,pvalue)
        names(result) = c("statistic","p.value")
        return(result)
    }
}

epsAC.rv.test.lmm = function(epsdata0,covariates0,RV,isxe,confounder,cind,weights){

    Bmat = diag(weights)

    y = epsdata0[,1]
    n = length(y)
    xg = as.matrix(RV)
    xe = as.matrix(covariates0)

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

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)
        egmat = matrix(eg,ncol=ng,nrow=n,byrow = TRUE)

        sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

        Sigma = sig1 + sig2
        s = (z[extreme]%*%g +sum(z[!extreme])*eg)/sigma2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)

        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        qdist = rep(0,1000000)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(1000000,1)
        }

        #critval = quantile(qdist,probs = 0.95)

        teststat = c(s%*%Bmat%*%t(s))

        pvalue = sum(qdist > teststat)/1000000

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
        coefs = fit$coef
        alpha = coefs[1]
        beta = coefs[2:length(coefs)]
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
        uindex = list()
        for(u in 1:nu){
            uindex[[u]] = (xe[extreme,cind] == ux[u,])
        }

        eg = list()
        varg =  list()
        egmat =  matrix(0,ncol=ng,nrow=n)

        for(u in 1:nu){
            uind = uindex[[u]]
            allind = (xe[,cind] == ux[u,])

            egmat[allind,] = colMeans(g[uind,])

            eg[[u]] = colMeans(g[uind,])
            varg[[u]] = var(g[uind,])
        }

        sig1 = (t(egmat)%*%(Ie-He)%*%egmat)/sigma2
        sig2 = 0

        s = z[extreme]%*%g

        for(u in 1:nu){
            uind = uindex[[u]]
            n_uc = length(z[extreme][uind])
            sig2 = sig2 + varg[[u]]*(sum(z[extreme][uind]^2) - (1/n_uc)*(sum(z[extreme][uind]))^2)/(sigma2^2)

            tempind = (xe[!extreme,cind] == ux[u,])
            s = s + sum(z[!extreme][tempind])*eg[[u]]
        }

        s = s/sigma2
        Sigma = sig1 + sig2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)
        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        qdist = rep(0,1000000)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(1000000,1)
        }

        #critval = quantile(qdist,probs = 0.95)

        teststat = c(s%*%Bmat%*%t(s))

        pvalue = sum(qdist > teststat)/1000000

        result = list(teststat,pvalue,d)
        names(result) = c("statistic","p.value","eigen")
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

        extreme = !is.na(xg[,1])
        g = as.matrix(xg[extreme,])
        n_cc = sum(extreme)
        ng = dim(g)[2]

        eg = colMeans(g)
        varg = var(g)

        sig2 = varg*(sum(z[extreme]^2) - (1/n_cc)*(sum(z[extreme]))^2)/(sigma2^2)

        Sigma = sig2
        s = (z[extreme]%*%g +sum(z[!extreme])*eg)/sigma2

        Vmat = eigen(Sigma)$vectors
        Dmat = t(Vmat)%*%Sigma%*%Vmat
        Smat = diag(sqrt(diag(Dmat)))
        sqrtSigma = Vmat%*%Smat%*%t(Vmat)
        eigenmat = sqrtSigma%*%Bmat%*%sqrtSigma
        d = eigen(eigenmat)$values

        qdist = rep(0,1000000)
        for(c in 1:length(d)){
            qdist = qdist + d[c]*rchisq(1000000,1)
        }

        #critval = quantile(qdist,probs = 0.95)

        teststat = c(s%*%Bmat%*%t(s))

        pvalue = sum(qdist > teststat)/1000000

        result = list(teststat,pvalue,d)
        names(result) = c("statistic","p.value","eigen")
        return(result)
    }
}
