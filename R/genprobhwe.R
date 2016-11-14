# returns all combinations of genotypes and their probability
# hardy-weinberg equilibrium is assumed

genprobhwe = function(ng,maf,geneffect = "additive"){
    temp1 = list()
    temp2 = list()
    if(geneffect == "additive"){
        for(i in 1:ng){
            temp1[[i]] = c(0,1,2)
            temp2[[i]] = c((1-maf[i])^2,2*maf[i]*(1-maf[i]),maf[i]^2)
        }
    }else if(geneffect == "dominant"){
        for(i in 1:ng){
            temp1[[i]] = c(0,1)
            temp2[[i]] = c((1-maf[i])^2,2*maf[i]*(1-maf[i])+maf[i]^2)
        }
    }else if(geneffect == "recessive"){
        for(i in 1:ng){
            temp1[[i]] = c(0,1)
            temp2[[i]] = c((1-maf[i])^2+2*maf[i]*(1-maf[i]),maf[i]^2)
        }
    }
    geno = expand.grid(temp1)
    probs = expand.grid(temp2)
    probs2 = apply(probs,1,prod)
    return(list(geno,probs2))
}
