## ------------------------------------------------------------------------
    library(extremesampling)

## ------------------------------------------------------------------------
    # Generate data set of the type EPS-only

    N = 5000 # Number of individuals in a population
    xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
    xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
    xg1 = sample(c(0,1,2),N,c(0.4,0.3,0.3), replace = TRUE) # SNP
    xg2 = sample(c(0,1,2),N,c(0.5,0.3,0.2), replace = TRUE) # SNP
    
    # Model parameters
    a = 50; be1 = 5; be2 = 8; bg1 = 0.3; bg2 = 0.6; sigma = 2
    
    # Generate response y (phenotype)
    y = rnorm(N, mean = a+be1*xe1+be2*xe2+bg1*xg1+bg2*xg2, sd = sigma)
    
    # Identify extremes, here upper and lower 25% of population
    u = quantile(y,probs = 3/4,na.rm=TRUE)
    l = quantile(y,probs = 1/4,na.rm=TRUE)
    extreme = (y < l) | (y >= u)
    
    # Create the EPS-only data set
    y = y[extreme]
    xe1 = xe1[extreme]
    xe2 = xe2[extreme]
    xg1 = xg1[extreme]
    xg2 = xg2[extreme]
    
    # Matrix form of covariates
    xg = as.matrix(cbind(xg1,xg2))
    xe = as.matrix(cbind(xe1,xe2))

## ------------------------------------------------------------------------
    # 1.1) Perform score test for one genetic covariate at a time
    epsonly.test(y~xe1+xe2, SNP=xg, cutoffs=c(l,u))$p.value
    # alternatively, all environmental covariates can be gathered in a matrix
    epsonly.test(y~xe, SNP=xg, cutoffs=c(l,u))$p.value
    
    # 1.2) Perform score test for all genetic covariates simultanously
    epsonly.test(y~xe1+xe2, SNP=xg, cutoffs=c(l,u), onebyone = FALSE)$p.value
    
    # 1.3) Test gene-environment interaction
    epsonly.testGE(y~xe1+xe2+xg1+xg2, GE=c("xe1:xg1"), cutoffs=c(l,u))$p.value

## ------------------------------------------------------------------------
    # 1.4) Fit a linear model
    epsonly.lm(y~xe1+xe2+xg1+xg2,cutoffs = c(l,u))
    # alternatively with matrix inputs
    # epsonly.lm(y~xe+xg,cutoffs = c(l,u))
    
    # 1.5) Fit a linear model with an interaction term
    epsonly.lm(y~xe1+xe2+xg1+xg2+xe2*xg1,cutoffs = c(l,u))$coefficients
    # alternatively with matrix inputs
    # epsonly.lm(y~xe+xg+xe2*xg1,cutoffs = c(l,u))$coefficients
    

## ------------------------------------------------------------------------
    # Generate example data set of the type EPS-complete

    N = 5000 # Number of individuals in a population
    xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
    xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
    xg1 = sample(c(0,1,2),N,c(0.4,0.3,0.3), replace = TRUE) # SNP
    xg2 = sample(c(0,1,2),N,c(0.5,0.3,0.2), replace = TRUE) # SNP
    
    # Model parameters
    a = 50; be1 = 5; be2 = 8; bg1 = 0.3; bg2 = 0.6; sigma = 2
    
    # Generate response y
    y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 + bg2*xg2, sd = sigma)
    
    # Identify extremes, here upper and lower 25% of population
    u = quantile(y,probs = 3/4,na.rm=TRUE)
    l = quantile(y,probs = 1/4,na.rm=TRUE)
    extreme = (y < l) | (y >= u)
    
    # Create the EPS-complete data set by setting
    # the SNP values of non-extremes to NA
    xg1[!extreme] = NA
    xg2[!extreme] = NA
    
    xg = as.matrix(cbind(xg1,xg2))
    xe = as.matrix(cbind(xe1,xe2))

## ------------------------------------------------------------------------
    # 2.1) Perform score test for one genetic covariate at a time
    epscomp.test(y~xe1+xe2,SNP = xg)$p.value
    # alternatively with matrix inputs
    epscomp.test(y~xe,SNP = xg)$p.value
    
    # 2.2) Perform score test for all genetic covariates simultanously
    epscomp.test(y~xe1+xe2,SNP = xg,onebyone = FALSE)$p.value
    
    # 2.3) Test interactions
    epscomp.testGE(y~xe1+xe2+xg1+xg2,GE = c("xe1:xg1"))$p.value
    # alternatively with matrix inputs
    # epscomp.testGE(y~xe+xg,GE = c("xe1:xg1"))$p.value

## ------------------------------------------------------------------------
    # 2.4) Fit a linear model
    epscomp.lm(y~xe1+xe2+xg1+xg2)
    # alternatively with matrix inputs
    # epscomp.lm(y~xe+xg)
    
    # 2.5) Fit a linear model with an interaction term
    epscomp.lm(y~xe1+xe2+xg1+xg2+xe1*xg2)
    # alternatively with matrix inputs
    # epscomp.lm(y~xe+xg+xe1*xg2)

## ------------------------------------------------------------------------
    # Generate example data set of the type EPS-complete
    # Assume Hardy-Weinberg

    N = 5000 # Number of individuals in a population
    xe = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
    q = 0.3 # minor allele frequency
    xg = sample(c(0,1,2),N,c((1-q)^2,2*q*(1-q),q^2), replace = TRUE) # SNP
    
    # Model parameters
    a = 50; be = 5; bg = 0.5; sigma = 2
    
    # Generate response y
    y = rnorm(N, mean = a + be*xe + bg*xg, sd = sigma)
    
    # Identify extremes, here upper and lower 25% of population
    u = quantile(y,probs = 3/4,na.rm=TRUE)
    l = quantile(y,probs = 1/4,na.rm=TRUE)
    extreme = (y < l) | (y >= u)
    
    # Create the EPS-complete data set by setting
    # the SNP values of non-extremes to NA
    xg[!extreme] = NA
    
    # Fit a linear model with known minor allele frequency
    epscomp.lm(y~xe+xg,hwe = TRUE, maf = q)
    
    # Fit a linear model with unknown minor allele frequency
    epscomp.lm(y~xe+xg,hwe = TRUE)

## ------------------------------------------------------------------------
    # Generate example data set of the type EPS-complete
    # with confounding

    N = 5000 # Number of individuals in a population
    xe1 = rnorm(n = N, mean = 2, sd = 1) # Environmental covariate
    xe2 = rbinom(n = N, size = 1, prob = 0.3) # Environmental covariate
    xg1 = sample(c(0,1,2),N,c(0.4,0.3,0.3), replace = TRUE) # SNP
    
    xg2 = rep(NA,N) # SNP
    n0 = sum(xe2 == 0)
    n1 = sum(xe2 == 1)
    xg2[xe2 == 0] = sample(c(0,1,2),n0,c(0.5,0.3,0.2), replace = TRUE)
    xg2[xe2 == 1] = sample(c(0,1,2),n1,c(0.3,0.1,0.6), replace = TRUE)
    
    # Model parameters
    a = 50; be1 = 5; be2 = 8; bg1 = 0.3; sigma = 2
    
    # Generate response y (xg2 has no effect)
    y = rnorm(N, mean = a + be1*xe1 + be2*xe2 + bg1*xg1 , sd = sigma)
    
    # Identify extremes, here upper and lower 25% of population
    u = quantile(y,probs = 3/4,na.rm=TRUE)
    l = quantile(y,probs = 1/4,na.rm=TRUE)
    extreme = (y < l) | (y >= u)
    
    xg1[!extreme] = NA
    xg2[!extreme] = NA
    
    xg = cbind(xg1,xg2)
    
    # Test (with confounding)
    epscomp.test(y~xe1+xe2,SNP = xg,confounder = TRUE, cx = c("xe2"))
    
    # Test (without confounding)
    epscomp.test(y~xe1+xe2,SNP = xg)

