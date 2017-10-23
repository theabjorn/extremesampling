# epsCC.test.W = function(nullmodel,SNP,y,yparam,randomindex){
#     options(na.action="na.pass")
#         epsdata0 = model.frame(nullmodel)
#         covariates0 = as.matrix(model.matrix(nullmodel)[,-1])
#     options(na.action="na.omit")
#
#     g = as.matrix(SNP)
#     ng = dim(g)[2]
#     totest = colnames(SNP)
#
#     statistic = matrix(NA,ncol = 1, nrow = ng)
#     pvalue = matrix(NA,ncol = 1, nrow = ng)
#     rownames(statistic) = totest
#     rownames(pvalue) = totest
#     colnames(statistic) = "t"
#     colnames(pvalue) = "p.value"
#
#     w = epsdata0[,1]
#     x = as.matrix(covariates0)
#     n = length(w)
#     gamma = (y - yparam[1] - cbind(x,g)%*%yparam[2:(length(yparam)-1)])/yparam[length(yparam)]
#
#     fit0 = epsCC.loglikmax.W(cbind(w,x),gamma,hessian = TRUE) # Fit under H0
#     fit = fit0[[2]]
#     alpha = fit[1]
#     beta = fit[2:(length(fit)-2)]
#     sigma = fit[length(fit)-1]
#     sigma2 = sigma*sigma
#     rho = fit[length(fit)]
#
#     xbeta = x%*%beta
#     z = w - alpha - xbeta - rho*sigma*gamma
#
#
#     I11_11 = (1/sigma2)*(1/(1-rho^2))*n
#     I11_22 = (1/sigma2)*(1/(1-rho^2))*t(x)%*%x
#     I11_21 = (1/sigma2)*(1/(1-rho^2))*colSums(x)
#     I11_12 = t(I11_21)
#
#     I11_31 = (1/sigma2)*(rho/(1-rho^2))*sum(gamma); I11_13 = t(I11_31)
#     I11_32 = (1/sigma2)*(rho/(1-rho^2))*t(gamma)%*%x; I11_23 = t(I11_32)
#
#     I11_41 = (1/sigma)*(1/(1-rho^2))*sum(gamma);  I11_14 = t(I11_41)
#     I11_42 = (1/sigma)*(1/(1-rho^2))*t(gamma)%*%x;  I11_24 = t(I11_42)
#
#     I11_33 = (1/sigma2)*(rho^2/(1-rho^2))*sum(gamma^2) - n/sigma^3
#     I11_44 = (1/(1-rho^2))*sum(gamma^2) - n*(1+rho^2)/((1-rho^2)^2)
#     I11_34 = (1/sigma2)*(rho/(1-rho^2))*sum(gamma^2)
#     I11_43 = t(I11_34)
#
#
#     I11 = cbind(c(I11_11,I11_21,I11_31,I11_41),
#                 rbind(I11_12,I11_22,I11_32,I11_42),
#                 rbind(I11_13,I11_23,I11_33,I11_43),
#                 rbind(I11_14,I11_24,I11_34,I11_44))
#
#     for(i in 1:ng){
#         gi = g[,i]
#
#         I22 = (1/sigma2)*(1/(1-rho^2))*t(gi)%*%gi
#
#         I21_1 = (1/sigma2)*(1/(1-rho^2))*sum(gi)
#         I21_2 = (1/sigma2)*(1/(1-rho^2))*t(gi)%*%x
#
#         I21_3 = (1/sigma2)*(1/(1-rho^2))*((2/sigma)*t(z)%*%gi + rho*t(gamma)%*%gi)
#         I21_4 = (1/sigma2)*(1/(1-rho^2))*(-2*rho/(1-rho^2)*t(z)%*%gi + sigma*t(gamma)%*%gi)
#
#         I21 = cbind(I21_1,I21_2,I21_3,I21_4)
#         I12 = t(I21)
#
#         Sigma = (I22 - I21%*%ginv(I11)%*%t(I21))
#         s = (1/sigma2)*(1/(1-rho^2))*t(gi)%*%(z)
#         t = s*s/Sigma
#         pval = pchisq(t,1,lower.tail=FALSE)
#
#         statistic[i,] = c(t)
#         pvalue[i,] = c(pval)
#     }
#     result = list(statistic,pvalue)
#     names(result) = c("statistic","p.value")
#     return(result)
# }
