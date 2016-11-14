# returns all combinations of genotypes and their probability

genprob = function(ng,probs,geneffect = "additive"){
    epsilon = 0.00001
    temp1 = list()
    temp2 = list()
    if(geneffect == "additive"){
        index = 1
        for(i in 1:ng){
            temp1[[i]] = c(0,1,2)
            temp2[[i]] = c(probs[index],probs[index+1],
                           (1-probs[index]-probs[index+1] + epsilon))
            index = index + 2
        }
    }else if(geneffect == "dominant"){
        for(i in 1:ng){
            temp1[[i]] = c(0,1)
            temp2[[i]] = c(probs[i],(1-probs[i] + epsilon))
        }
    }else if(geneffect == "recessive"){
        for(i in 1:ng){
            temp1[[i]] = c(0,1)
            temp2[[i]] = c(probs[i],(1-probs[i]))
        }
    }
    geno = expand.grid(temp1)
    probs = expand.grid(temp2)
    probs2 = apply(probs,1,prod)
    return(list(geno,probs2))
}
