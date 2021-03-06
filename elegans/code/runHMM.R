# emissionProbs is a list of sparse matrices with emission probabilities given observed count data for parental variants (pRR, pHet, and pAA)
# hmm.out.dir , where to save info for each chromosome
# gmap.s, genetic map, pre-jittered so that genetic distance between 

runHMM=function(emissionProbs, hmm.out.dir, gmap.s, chroms,  calc='genoprob',n.indiv=NULL, sex.chr.model=F,lik=F) {
    attach(emissionProbs)
    #pRR
    #pHet
    #pAA
    #calc='genoprobs'
    #calc='viterbi'
    #n.indiv for debugging
    if(is.null(n.indiv)) {    n.indiv=ncol(pRR) }
    sex.chr=F
    for(coii in chroms) {
        if(coii=='X' & sex.chr.model) {sex.chr=T}

        cell.names=colnames(pRR)
       
        coi=paste0('^', coii, '_')
        # marker indices for relevant chromosome 
        mind=grep(coi, rownames(pRR))

        #genetic map
        gdist=gmap.s[[coii]]$map 
        names(gdist)=gmap.s[[coii]]$marker 
        
        #convert to recombination fraction
        g.rf=map2rf(gdist)
        #convert to transmission matrices 
        tMats=lapply(g.rf, mTmat)

        #double check this 
        #if(sex.chr){ tMats=lapply(g.rf, mTmatH) }
        if(sex.chr){ startProbs=c(.44,.44,.11) } else {startProbs=NULL } 

        #initialize storage for posterior probabilities
        split.indiv=cut(1:n.indiv, ncores) 
        
        #relying on .combine is terribly slow, so don't
        pp=foreach(subset.indiv=levels(split.indiv)) %dopar% { 
            indiv.oi=which(split.indiv==subset.indiv)
            if(calc=='genoprob') {
                posteriorProb=array(0, dim=c(length(indiv.oi),3,length(gdist)),
                                   dimnames=list(id=cell.names[indiv.oi],
                                                 geno=c('AA','AB', 'BB'),
                                                 markers=names(gdist)) )
                #if(sex.chr) {posteriorProb=posteriorProb[,-2,] }
                lik.vec=rep(NA,length(indiv.oi))
                names(lik.vec)=cell.names[indiv.oi]
            }  
            if(calc=='viterbi') {
                    viterbiOut=matrix(0, length(indiv.oi),  length(gdist),
                                  dimnames=list(id=cell.names[indiv.oi],
                                                markers=names(gdist)) )                        
            }
            for(indiv in indiv.oi) { 
                if(indiv%%100==0) {print(paste(coii, indiv))}
                indRR=pRR[mind,indiv]
                indHet=pHet[mind,indiv]
                indAA=pAA[mind,indiv]
                eMats=rbind(indRR,indHet,indAA)
                ess=colSums(eMats)
                eMats[,ess<1e-12]=1
                #if(sex.chr) {eMats=eMats[-2,]}

                if(calc=='genoprob') {
                    #print(startProbs)
                    f=calcForward(eMats,tMats,startProbs)
                    b=calcBackward(eMats,tMats)
                    posteriorProb[cell.names[indiv],,]=calcPosterior(f,b)
                    if(lik) { lik.vec[cell.names[indiv]]=calcForward(posteriorProb[cell.names[indiv],,],tMats,lik=T)}
                }
                if(calc=='viterbi') {
                    print(startProbs)
                    test=viterbi(eMats,tMats,startProbs)
                    viterbiOut[cell.names[indiv],]=viterbi(eMats,tMats)
                 }
            }
            if(lik){return(lik.vec)}
            if(calc=='genoprob') {  return(posteriorProb) }
            if(calc=='viterbi')   {  return(viterbiOut)    } 
        }
       # het.liks=do.call('c', pp)
       # hap.liks=do.call('c', pp)
       # pchisq(-2*((hap.liks)-(het.liks)),1, lower.tail=F)
       # plot(hap.liks,het.liks, col=ifelse(pchisq(-2*((hap.liks)-(het.liks)),1, lower.tail=F)<.05,'red','black'))
        
        
        if(calc=='genoprob') { 
            #faster to do it this way then rely on foreach .combine  
            posteriorProb=array(0, dim=c(n.indiv,3,length(gdist)),
                                   dimnames=list(id=cell.names[1:n.indiv],
                                                 geno=c('AA','AB', 'BB'),
                                                 markers=names(gdist))   )
            #if(sex.chr) {posteriorProb=posteriorProb[,-2,] }

            # and then fill it up manually with a loop
            for(n in 1:length(pp)) {
                #print(n)
                posteriorProb[rownames(pp[[n]]),,]=pp[[n]]
            }
            saveRDS(posteriorProb, file=paste0(hmm.out.dir, coii,'.RDS'))
        }
        if(calc=='viterbi') {
            viterbiOut=matrix(0, n.indiv,  length(gdist),
                                  dimnames=list(id=cell.names[1:n.indiv],
                                                markers=names(gdist)) ) 
         for(n in 1:length(pp)) {
                print(n)
                viterbiOut[rownames(pp[[n]]),]=pp[[n]]
            }

            saveRDS(viterbiOut, file=paste0(hmm.out.dir, coii,'_viterbi.RDS'))
        }
       # return(NULL)
    }
    return(NULL)
}
