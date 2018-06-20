# if(hwe){
#     #######################################################################
#     # 1: Hardy-Weinberg
#     #######################################################################
#     if(interact == FALSE){
#         ###################################################################
#         # 1.1: No interactions
#         ###################################################################
#         if(!isxe){
#             ###############################################################
#             # 1.1.1: No xe
#             ###############################################################
#             if(!is.na(maf)){
#                 ###########################################################
#                 # 1.1.1.1: MAF given
#                 ###########################################################
#                 message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequency: ", toString(maf),sep = ""))
#                 data = cbind(y,xg)
#                 ng = dim(as.matrix(xg))[2]
#                 model = epsAC.loglikmax(data, ng, hwe = TRUE, maf = maf,
#                                         hessian = TRUE)
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 params = model[[2]]
#                 sigma = params[length(params)]
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = (length(params)-1),
#                                        ncol = 2))
#                 for(i in 1:(length(params)-1)){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)", newnames)
#                 rownames(ci) = c("(intercept)",newnames)
#                 result = list(coef,ci,sigma)
#                 names(result) = c("coefficients","ci","sigma")
#                 return(result)
#             }else{
#                 ###########################################################
#                 # 1.1.1.2: MAF not given
#                 ###########################################################
#                 message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
#                 data = cbind(y,xg)
#                 ng = dim(as.matrix(xg))[2]
#                 model = epsAC.loglikmax(data, ng, hwe = TRUE,
#                                         hessian = TRUE)
#                 params = model[[2]]
#                 nparam = 1 + ng
#                 sigma = params[(nparam + 1)]
#                 maf = params[(nparam + 2):length(params)]
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = (length(params)-1-ng),
#                                        ncol = 2))
#                 for(i in 1:nparam){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)",newnames)
#                 rownames(ci) = c("(intercept)",newnames)
#                 result = list(coef,ci,sigma,maf)
#                 names(result) = c("coefficients","ci","sigma","maf")
#                 return(result)
#             }
#         }else{
#             ###############################################################
#             # 1.1.2: xe
#             ###############################################################
#             data = cbind(y,xe,xg)
#             xe = as.matrix(xe)
#             xg = as.matrix(xg)
#             ng = dim(as.matrix(xg))[2]
#             ne = dim(as.matrix(xe))[2]
#
#             if(!is.na(as.matrix(maf)[1,1])){
#                 ###########################################################
#                 # 1.1.2.1: MAF given
#                 ###########################################################
#                 message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequencies: ", toString(maf),sep = ""))
#                 model = epsAC.loglikmax(data,ng, hwe = TRUE, maf = maf,
#                                         hessian = TRUE)
#                 params = model[[2]]
#                 nparam = 1 + ne + ng
#                 sigma = params[(nparam + 1)]
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = (1 + ne + ng), ncol = 2))
#                 for(i in 1:(1 + ne + ng)){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)",newnames)
#                 rownames(ci) = c("(intercept)",newnames)
#
#                 result = list(coef,ci,sigma)
#                 names(result) = c("coefficients","ci","sigma")
#                 return(result)
#             }else{
#                 ###########################################################
#                 # 1.1.2.2: MAF not given
#                 ###########################################################
#                 message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
#                 model = epsAC.loglikmax(data, ng, hwe = TRUE,
#                                         hessian = TRUE)
#                 params = model[[2]]
#                 nparam = 1 + ne + ng
#                 sigma = params[(nparam + 1)]
#                 maf = params[(nparam + 2):length(params)]
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = (length(params)-1-ng),
#                                        ncol = 2))
#                 for(i in 1:(1 + ne + ng)){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)",newnames)
#                 rownames(ci) = c("(intercept)",newnames)
#                 result = list(coef,ci,sigma,maf)
#                 names(result) = c("coefficients","ci","sigma","maf")
#                 return(result)
#             }
#         }
#     }else{
#         ###################################################################
#         # 1.2: Interactions
#         ###################################################################
#         data = cbind(y,xe,xg)
#         xe = as.matrix(xe)
#         xg = as.matrix(xg)
#         ng = dim(xg)[2]
#         ne = dim(xe)[2]
#         intnames = c()
#         for(i in 1:length(interactind)){
#             intnames[i] = paste(colnames(xe)[interactind[[i]][2]],"*",
#                                 colnames(xg)[interactind[[i]][1]],sep="")
#         }
#         if(!is.na(as.matrix(maf)[1,1])){
#             ###########################################################
#             # 1.2.1: MAF given
#             ###########################################################
#             message(paste("Hardy-Weinberg equilibrium assumed with known minor allele frequency: ", toString(maf),sep = ""))
#             model = epsAC.loglikmaxint(data, ng, interactind,
#                                        hwe = TRUE, maf = maf,
#                                        hessian = TRUE)
#             params = model[[2]]
#             nparam = 1 + ne + ng + length(interactind)
#             sigma = params[(nparam + 1)]
#             hessian = model[[1]]
#             info = -1*ginv(hessian)
#             coef = c()
#             ci = data.frame(matrix(NA,nrow = (nparam), ncol = 2))
#             for(i in 1:(nparam)){
#                 coef[i] = params[i]
#                 ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                 ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#             }
#             colnames(ci) = c("lower 95% ci", "upper 95% ci")
#             names(coef) = c("(intercept)",newnames)
#             rownames(ci) = c("(intercept)",newnames)
#             result = list(coef,ci,sigma)
#             names(result) = c("coefficients","ci","sigma")
#             return(result)
#         }else{
#             ###########################################################
#             # 1.2.1: MAF not given
#             ###########################################################
#             message("Hardy-Weinberg equilibrium assumed, unknown minor allele frequency.")
#             model = epsAC.loglikmaxint(data, ng, interactind,
#                                        hwe = TRUE,
#                                        hessian = TRUE)
#             params = model[[2]]
#             nparam = 1 + ne + ng + length(interactind)
#             sigma = params[(nparam + 1)]
#             maf = params[(nparam + 2):length(params)]
#             hessian = model[[1]]
#             info = -1*ginv(hessian)
#             coef = c()
#             ci = data.frame(matrix(NA,nrow = (nparam), ncol = 2))
#             for(i in 1:(nparam)){
#                 coef[i] = params[i]
#                 ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                 ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#             }
#             colnames(ci) = c("lower 95% ci", "upper 95% ci")
#             names(coef) = c("(intercept)",newnames)
#             rownames(ci) = c("(intercept)",newnames)
#             result = list(coef,ci,sigma,maf)
#             names(result) = c("coefficients","ci","sigma","maf")
#             return(result)
#         }
#     }
# }else{
#     #####################################################################
#     # 2. Hardy-Weinberg not assumed
#     #####################################################################
#     if(interact == FALSE){
#         #################################################################
#         # 2.1.1: No interactions
#         #################################################################
#         if(!isxe){
#             ###############################################################
#             # 2.1.1: No xe
#             ###############################################################
#             data = cbind(y,xg)
#             ng = dim(as.matrix(xg))[2]
#             model = epsAC.loglikmax(data,ng,hessian = TRUE)
#             params = model[[2]]
#             nparam = 1 + ng
#             sigma = params[(nparam + 1)]
#             hessian = model[[1]]
#             info = -1*ginv(hessian)
#             coef = c()
#             ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
#             for(i in 1:nparam){
#                 coef[i] = params[i]
#                 ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                 ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#             }
#             gfreqs = model[[3]]
#             genotypes = model[[4]]
#             resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
#             colnames(resg) = c("P(Xg)",covnames[snpid])
#             colnames(ci) = c("lower 95% ci", "upper 95% ci")
#             names(coef) = c("(intercept)",newnames)
#             rownames(ci) = c("(intercept)",newnames)
#             result = list(coef,ci,sigma,resg)
#             names(result) = c("coefficients","ci","sigma","Xg")
#             return(result)
#         }else{
#             ###############################################################
#             # 2.1.1: xe
#             ###############################################################
#             data = cbind(y,xe,xg)
#             xe = as.matrix(xe)
#             xg = as.matrix(xg)
#             ng = dim(as.matrix(xg))[2]
#             ne = dim(as.matrix(xe))[2]
#             if(!confounder){
#                 ###########################################################
#                 # 2.1.1: No confounders
#                 ###########################################################
#                 model = epsAC.loglikmax(data, ng,hessian = TRUE)
#                 params = model[[2]]
#                 nparam = 1 + ne + ng # alpha, betae, betag
#                 sigma = params[(nparam + 1)]
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
#                 for(i in 1:(1 + ne + ng)){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 gfreqs = model[[3]]
#                 genotypes = model[[4]]
#                 resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
#                 colnames(resg) = c("P(Xg)",covnames[snpid])
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)",newnames)
#                 rownames(ci) = c("(intercept)",newnames)
#                 result = list(coef,ci,sigma,resg)
#                 names(result) = c("coefficients","ci","sigma","Xg")
#                 return(result)
#             }else{
#                 #########################################################
#                 # 2.1.1: Confounders
#                 #########################################################
#                 model = epsAC.loglikmaxcond(data, ng,cind = cind,
#                                             hessian = TRUE,snpnames=covnames[snpid])
#                 params = model[[2]]
#                 nparam = 1 + ne + ng
#                 sigma = params[(nparam + 1)]
#                 probs = params[(nparam + 2):length(params)]
#                 hessian = model[[1]]
#                 info = -1*ginv(hessian)
#                 coef = c()
#                 ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
#                 for(i in 1:(1 + ne + ng)){
#                     coef[i] = params[i]
#                     ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                     ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#                 }
#                 resg = model[[3]]
#                 colnames(ci) = c("lower 95% ci", "upper 95% ci")
#                 names(coef) = c("(intercept)",
#                                 newnames)
#                 rownames(ci) = c("(intercept)",
#                                  newnames)
#                 result = list(coef,ci,sigma,resg)
#                 names(result) = c("coefficients","ci","sigma","Xg")
#                 return(result)
#             }
#         }
#     }else{
#         #################################################################
#         # 2.1.1: Interactions
#         #################################################################
#         data = cbind(y,xe,xg)
#         xe = as.matrix(xe)
#         xg = as.matrix(xg)
#         ng = dim(as.matrix(xg))[2]
#         ne = dim(as.matrix(xe))[2]
#         intnames = c()
#         for(i in 1:length(interactind)){
#             intnames[i] = paste(colnames(xe)[interactind[[i]][2]],"*",
#                                 colnames(xg)[interactind[[i]][1]],sep="")
#         }
#         if(!confounder){
#             ###############################################################
#             # 2.1.1: No confounders
#             ###############################################################
#             model = epsAC.loglikmaxint(data, ng, interactind = interactind,hessian = TRUE)
#             params = model[[2]]
#             nparam = 1 + ne + ng + length(interactind)
#             sigma = params[(nparam + 1)]
#             hessian = model[[1]]
#             info = -1*ginv(hessian)
#             coef = c()
#             ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
#             for(i in 1:nparam){
#                 coef[i] = params[i]
#                 ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                 ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#             }
#             gfreqs = model[[3]]
#             genotypes = model[[4]]
#             resg = data.frame(cbind(gfreqs,genotypes),row.names = NULL)
#             colnames(resg) = c("P(Xg)",covnames[snpid])
#             colnames(ci) = c("lower 95% ci", "upper 95% ci")
#             names(coef) = c("(intercept)",newnames)
#             rownames(ci) = c("(intercept)",newnames)
#             result = list(coef,ci,sigma,resg)
#             names(result) = c("coefficients","ci","sigma","Xg")
#             return(result)
#         }else{
#             ###############################################################
#             # 2.1.1: Confounders
#             ###############################################################
#             model = epsAC.loglikmaxcondint(data, ng, cind = cind,
#                                            interactind = interactind,
#                                            hessian = TRUE,snpnames=covnames[snpid])
#             params = model[[2]]
#             nparam = 1 + ne + ng + length(interactind)
#             sigma = params[(nparam + 1)]
#             hessian = model[[1]]
#             info = -1*ginv(hessian)
#             coef = c()
#             ci = data.frame(matrix(NA,nrow = nparam, ncol = 2))
#             for(i in 1:nparam){
#                 coef[i] = params[i]
#                 ci[i,1] = params[i] - 1.96*(sqrt(info[i,i]))
#                 ci[i,2] = params[i] + 1.96*(sqrt(info[i,i]))
#             }
#             resg = model[[3]]
#             colnames(ci) = c("lower 95% ci", "upper 95% ci")
#             names(coef) = c("(intercept)",
#                             newnames)
#             rownames(ci) = c("(intercept)",
#                              newnames)
#             result = list(coef,ci,sigma,resg)
#             names(result) = c("coefficients","ci","sigma","Xg")
#             return(result)
#         }
#     }
# }
