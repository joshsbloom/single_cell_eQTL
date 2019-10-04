 estimateEmissionProbs=function(N2.counts, CB.counts, error.rate=.002) {
    #generates genotype emission probabilities as pRR, pHet, and pAlt
    # use error model from gusmap to estimate genotype emission probabilities
    # 10.1534/genetics.117.300627
    # eq(12) and SuppFile 1

    n=N2.counts+CB.counts
    k=N2.counts
    # per UMI genotype-specific count error rate
    ee=error.rate

    # simpler to keep this as a vector for these functions
    hascounts=which(n>0)
    pRR  = dbinom(n[hascounts]-k[hascounts],n[hascounts],ee) 
    pHet = choose(n[hascounts],k[hascounts])*(1/(2^n[hascounts]))
    pAA  = dbinom(k[hascounts],n[hascounts],ee) 

    # simpler to keep this as an array for sparseMatrix()
    hascounts=which(n>0,arr.ind=T) 
    pRR  = sparseMatrix(i=hascounts[,1], j=hascounts[,2], x=pRR,  dims=dim(n), dimnames=dimnames(n))
    pHet = sparseMatrix(i=hascounts[,1], j=hascounts[,2], x=pHet, dims=dim(n), dimnames=dimnames(n))
    pAA  = sparseMatrix(i=hascounts[,1], j=hascounts[,2], x=pAA,  dims=dim(n), dimnames=dimnames(n))

    # check these bastards, a couple hundred instances of NAs coming from unreliable or super high N sites
    pHet[which(is.na(pHet),arr.ind=T)]=0
    return(list(pRR=pRR,pHet=pHet,pAA=pAA))
    
    # fix code to make downstream code not depend on this !
    # especially, I_5019104 indiv 14,20,174
    # nh=n[hascounts]
    # kh=k[hascounts]
}


