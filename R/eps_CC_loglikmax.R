# Maximize log-likelihood function for EPS-only data

epsCC.loglikmax = function(data,l,u,randomindex,ll=FALSE, hessian = FALSE){

    rsample = TRUE
    if(missing(randomindex)){
        rsample = FALSE
    }else if(sum(randomindex)==0){
        rsample = FALSE
    }

    if(!rsample){
        if(dim(data)[2]>1){
            y = data[,1]
            len = dim(data)[2]
            x = as.matrix(data[,2:len])
            init = lm(y ~ x)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik_ex,
                           data = data,
                           l=l,
                           u=u,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }else{
            y = data[,1]
            init = lm(y ~ 1)
            a.init = c(init$coef,log(summary(init)$sigma^2))
            result = optim(a.init,
                           epsCC.loglik_e,
                           data = data,
                           l=l,
                           u=u,
                           method = "BFGS",
                           control = list(fnscale = -1,
                                          trace = 1),
                           hessian = TRUE)
            if(ll){
                return(result$value)
            }else if(hessian){
                return(list(result$hessian,
                            c(result$par[1:(length(result$par)-1)],
                              sqrt(exp(result$par[(length(result$par))])))))
            }else{
                return(c(result$par[1:(length(result$par)-1)],
                         sqrt(exp(result$par[(length(result$par))]))))
            }
        }

    }else{
        # if(dim(data)[2]>1){
        #     y = data[,1]
        #     len = dim(data)[2]
        #     x = as.matrix(data[,2:len])
        #     init = lm(y ~ x)
        #     a.init = c(init$coef,log(summary(init)$sigma^2))
        #     result = optim(a.init,
        #                    epsCC.loglik_rand,
        #                    data = data,
        #                    l=l,
        #                    u=u,
        #                    randomindex = randomindex,
        #                    method = "BFGS",
        #                    control = list(fnscale = -1,
        #                                   trace = 1),
        #                    hessian = TRUE)
        #     if(ll){
        #         return(result$value)
        #     }else if(hessian){
        #         return(list(result$hessian,
        #                     c(result$par[1:(length(result$par)-1)],
        #                       sqrt(exp(result$par[(length(result$par))])))))
        #     }else{
        #         return(c(result$par[1:(length(result$par)-1)],
        #                  sqrt(exp(result$par[(length(result$par))]))))
        #     }
        # }else{
        #     y = data[,1]
        #     init = lm(y ~ 1)
        #     a.init = c(init$coef,log(summary(init)$sigma^2))
        #     result = optim(a.init,
        #                    epsCC.loglik_z_rand,
        #                    data = data,
        #                    l=l,
        #                    u=u,
        #                    randomindex=randomindex,
        #                    method = "BFGS",
        #                    control = list(fnscale = -1,
        #                                   trace = 1),
        #                    hessian = TRUE)
        #     if(ll){
        #         return(result$value)
        #     }else if(hessian){
        #         return(list(result$hessian,
        #                     c(result$par[1:(length(result$par)-1)],
        #                       sqrt(exp(result$par[(length(result$par))])))))
        #     }else{
        #         return(c(result$par[1:(length(result$par)-1)],
        #                  sqrt(exp(result$par[(length(result$par))]))))
        #     }
        # }
    }
}




# #######################################################
# # EPS-CC loglik for secondary phenotype W
# #######################################################
#
# epsCC.loglikmax.W = function(data,gamma,ll=FALSE, hessian = FALSE){
#     len = dim(data)[2]
#     w = data[,1]
#     x = as.matrix(data[,2:(len)])
#     init = lm(w ~ x)
#     rho.init = 0.5
#     a.init = c(init$coef,log(summary(init)$sigma^2),log(rho.init/(1-rho.init)))
#
#     data = cbind(w,x,gamma)
#     result = optim(a.init,
#                    epsCC.loglik.W,
#                    data = data,
#                    method = "BFGS",
#                    control = list(fnscale = -1,
#                                   trace = 1),
#                    hessian = TRUE)
#     if(ll){
#         return(result$value)
#     }else if(hessian){
#         return(list(result$hessian,
#                     c(result$par[1:(length(result$par)-2)],
#                       sqrt(exp(result$par[(length(result$par-1))])),
#                       exp(result$par[(length(result$par))])/(1+exp(result$par[(length(result$par))])))))
#     }else{
#         return(c(result$par[1:(length(result$par)-2)],
#                  sqrt(exp(result$par[(length(result$par-1))])),
#                  exp(result$par[(length(result$par))])/(1+exp(result$par[(length(result$par))]))))
#     }
# }
