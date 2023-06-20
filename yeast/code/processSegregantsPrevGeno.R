 getSegMatch=function(comb.out.dir, vg, m.granges, geno_lookup_file='/data/eQTL/gdata_42k.RData') {
 

        load(geno_lookup_file) #'/data/eQTL/gdata_42k.RData')
        prevG=gdata
        rm(gdata)
        colnames(prevG)=gsub(':', '_', colnames(prevG))
        scname=tstrsplit(colnames(prevG),'_')
        scname=paste0(scname[[1]], '_', scname[[2]], '_', tstrsplit(scname[[3]],'/')[[1]], '_',tstrsplit(scname[[3]],'/')[[2]] )
        colnames(prevG)=scname

    # line up with existing genotypes 
        vg2=vg
        vg2name=tstrsplit(colnames(vg2), '_')
        colnames(vg2)=paste0(vg2name[[1]], '_', vg2name[[2]], '_', vg2name[[3]], '_', vg2name[[4]])

        matchingmarkers=colnames(vg2)[colnames(vg2) %in% colnames(prevG)]
        vg2=vg2[,matchingmarkers]

        pruned2=LDprune(vg2, m.granges)
        markerGRr2=pruned2$markerGRr
      
        prevGr=prevG[,matchingmarkers]
        stvg2=scale(t(vg2))
        sprevG=scale(t(prevGr))   
 
        # find best match genotype
        allcors=crossprod(stvg2, sprevG)/(nrow(sprevG)-1)
        best_match_seg=rownames(prevG)[(apply(allcors, 1 , which.max))]
        best_match_seg.r=apply(allcors, 1 , max)
 
        segMatch=data.frame(mmp1, best_match_seg=best_match_seg, best_match_seg.r=best_match_seg.r)
        saveRDS(segMatch, file=paste0(comb.out.dir, 'segMatch.RDS'))

        gdataPrev=(prevG+1)/2

        return(list(best_match_seg=best_match_seg,
                    gdataPrev=gdataPrev,
                    sprevG=sprevG,
                    markerGRr2=markerGRr2
                    ))                    
 
 }
 dobGLMM=function(comb.out.dir, data, Yr) {
        bGLMMs=list()
 #check out 2753
 #check out 3121
 #check out 3724
 #check out 4110
        for(i in 1:ncol(Yr)){
           print(i)
           tryCatch({
           bGLMMs[[colnames(Yr)[i]]]=glmmTMB(Yr[,i]~l_ctot+expt+(1|Cid)+(1|Zid), family=nbinom2, data=data, control=glmmTMBControl(parallel=36))
           },error=function(e) { next;
       }
       saveRDS(bGLMMs, file=paste0(comb.out.dir, 'bGLMMs.RDS'))
}        

dobGLMMtesting=function(comb.out.dir, data, Yr) {
        bGLMMs=list()
  
         for(i in 1:ncol(Yr)){
            print(i)
            bGLMMs[[colnames(Yr)[i]]]=glmmTMB(Yr[,i]~l_ctot+expt+(1|Cid)+(1|Zid), family=nbinom2, data=data, control=glmmTMBControl(parallel=36))
        }
        saveRDS(bGLMMs, file=paste0(comb.out.dir, 'bGLMMs.RDS'))
        

    sus=sort(unique(best_match_seg))
    A=tcrossprod(scale(t(sprevG[,sus])))/nrow(sprevG)
    A2=data.matrix(nearPD(A,corr=T)$mat)
    A2.inv=solve(A2)

    data2=data

    bGLMMsA=list()
    for( i in sigY ){
    #f= #Poisson(link="log"))
        print(colnames(Yr)[i])
        bGLMMsA[[colnames(Yr)[i]]]=

            system.time({aa=spaMM::fitme(Yr[,i]~l_ctot+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='REML', data=data2, family=negbin2)})
            system.time({bb=spaMM::fitme(Yr[,i]~l_ctot+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2)})
            system.time({bb=spaMM::fitme(Yr[,i]~l_ctot+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), data=data2, family=negbin2)})
            
            system.time({cc=spaMM::fitme(Yr[,i]~offset(l_ctot)+expt+(1|Zid)+(1|Cid), method='REML', data=data2, family=negbin2)})
            system.time({ff=spaMM::fitme(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), method='PQL', data=data2, family=negbin2)})
            system.time({ff=spaMM::fitme(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), method='REML', data=data2, family=negbin2, nb_cores=72)})
            system.time({cc=spaMM::corrHLfit(Yr[,i]~l_ctot+expt+corrFamily(1|Zid)+corrFamily(1|Cid), data=data2, family=negbin2)})

            system.time({cc=spaMM::fitme(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), method='ML',  data=data2, family=negbin2)})
            system.time({dd=spaMM::fitme(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), method='REML', data=data2, family=negbin2)})
            system.time({ee=glmmTMB(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), family=nbinom2, data=data2, REML=T, control=glmmTMBControl(parallel=36))})
            system.time({gg=vglmer(Yr[,i]~l_ctot+expt+(1|Zid)+(1|Cid), family=negbin, data=data2)}) #, REML=T, control=glmmTMBControl(parallel=36))})
            
            VarCorr(ee)$cond$Zid[1]
            VarCorr(ee)$cond$Cid[1]


    
            print(bGLMMsA[[colnames(Yr)[i]]])


    }

 }


 doLocalTestPrevGeno=function(sgd.genes, markerGRr2, Yr, dispersion.df, Gsub, gdataPrev, comb.out.dir, cl) {
        cisMarkers2=getCisMarker(sgd.genes[sgd.genes$gene_id %in% dispersion.df$gene,],
                            markerGRr2, t(Yr[,colnames(Yr) %in% dispersion.df$gene]) ) 
        n=tstrsplit(cisMarkers2$marker, '_')[-5]

        cisMarkers2$marker=paste0(n[[1]],'_', n[[2]], '_', n[[3]], '_', n[[4]])
        cisMarkers2=cisMarkers2[-which(is.na(match(cisMarkers2$marker, colnames(gdataPrev)))),]
        
        #counts)
        cSplit2=split(cisMarkers2, cisMarkers2$marker)

        GsubPlaceholder=Gsub
        cSplitPlacheholder=cSplit

        Gsub=gdataPrev[best_match_seg, ] #as.character(df.id),]
        rownames(Gsub)=rownames(GsubPlaceholder) #[order(df$id)]
        cSplit=cSplit2
        Gsub=Gsub[, unique(cisMarkers2$marker)]

        clusterExport(cl, varlist=c("mmp1", "Gsub", "counts","best_match_seg","doCisNebula3"))  #, "nebula"))
        system.time({   cisNB=parLapply(cl, cSplit, doCisNebula3)})
        names(cisNB)=names(cSplit)
        
        saveRDS(cisNB, file=paste0(comb.out.dir, 'cisNB2_prevG.RDS'))
        cis.ps=as.vector(unlist(sapply(cisNB, function(x) x$summary$p_)))
        cis.beta=as.vector(unlist(sapply(cisNB, function(x) x$summary$logFC_)))
        cis.gene=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))

        pGenoCis=data.frame(gene=cis.gene, beta=cis.beta, p=cis.ps)
        return(pGenoCis)
 }

calcICCPrevSegs=function(comb.out.dir, Yr) {

    bGLMMs=readRDS(paste0(comb.out.dir, 'bGLMMs.RDS'))

    #llikF=rep(0,length(bGLMMs))
    llikF=sapply(bGLMMs, logLik)
    llikmZ=rep(NA,length(bGLMMs))
    llikmC=rep(NA, length(bGLMMs))
    names(llikmZ)=names(llikF)
    names(llikmC)=names(llikF)

    #fit reduced model 
    for(i in 1:ncol(Yr)){
        if(is.na(llikF[i])){next;}
        print(i)
        x=bGLMMs[[i]]
        llikmZ[i]= logLik(update(x, . ~.-(1|Zid),control=glmmTMBControl(parallel=36)))
        llikmC[i]= logLik(update(x, . ~.-(1|Cid),control=glmmTMBControl(parallel=36)))
    }
    scVCmodels=data.frame(gene=names(llikF),llikF=llikF, llikmZ=llikmZ, llikmC=llikmC)
    saveRDS(scVCmodels, file=paste0(comb.out.dir, 'scVCmodels.RDS'))
    
    scVCmodels=readRDS(paste0(comb.out.dir, 'scVCmodels.RDS'))

    #table of ICC results #check out 2753
 #check out 3121
 #check out 3724
 #check out 4110

    ssH2=matrix(NA,length(bGLMMs),10)
    for( i in 1:length(bGLMMs) ) { #ncol(Yr)) { #[-c(2753,3121,3724,4110)]) {    
    #x=bGLMMs[[500]]
    #performance::icc(x)
        print(i)
        x=bGLMMs[[i]]
        lmbda=mean(exp((predict(x))))
        #print(lmbda)
        sgma=sigma(x)
        ssH2[i,8]=sgma
        ssH2[i,9]=lmbda
        ssH2[i,10]=mean(x$frame[[1]])
        sA=as.numeric(VarCorr(x)$cond$Zid[1]) #VarCorr(x)$cond$Zid[1]
        sC=as.numeric(VarCorr(x)$cond$Cid[1]) #VarCorr(x)$cond$Cid[1]
        ssH2[i,1]=sC/(sA+sC+log(1+(1/lmbda)+(1/sgma)))
        ssH2[i,2]=sA/(sA+sC+log(1+(1/lmbda)+(1/sgma)))
        ssH2[i,3]=sC
        ssH2[i,4]=sA

     
    #    trigamma(1/((1/lmbda)+(1/sgma)))
    }
    rownames(ssH2)=names(bGLMMs) #colnames(Yr)
    #ssH2[is.na(llikF),]=NA

    ssH2[,5]=p.adjust(pchisq(-2*(scVCmodels[names(bGLMMs),]$llikmC-scVCmodels[names(bGLMMs),]$llikF),1,lower.tail=F),method='fdr')
    ssH2[,6]=p.adjust(pchisq(-2*(scVCmodels[names(bGLMMs),]$llikmZ-scVCmodels[names(bGLMMs),]$llikF),1, lower.tail=F), method='fdr')
    ssH2[,7]=tot.umis=apply(Yr[,names(bGLMMs)],2,sum,na.rm=T)
    colnames(ssH2)=c('ICC.cc', 'ICC.H2', 'ccVar', 'H2Var', 'CC.q', 'H2.q', 'tot.umis', 'sigma', 'lambda', 'mean_obs')
    saveRDS(ssH2, file=paste0(comb.out.dir, 'scICC_2.RDS'))
    return(ssH2)
}


#take subset of genes, take previous lookup genotypes, and fit narrow-sense heritability
calcICCPrevSegs_3VC=function(comb.out.dir , Yr,sigY, best_match_seg, sprevG, data) {
       
    #comb.out.dir='/data/single_cell_eQTL/yeast/results/combined/Ap/'               
    
    sus=sort(unique(best_match_seg))
    A=tcrossprod(scale(t(sprevG[,sus])))/nrow(sprevG)
    A2=data.matrix(nearPD(A,corr=T)$mat)
    A2.inv=solve(A2)

    data2=data

    bGLMMsA=list()
    for( i in sigY ){
    #f= #Poisson(link="log"))
        print(colnames(Yr)[i])
        bGLMMsA[[colnames(Yr)[i]]]=spaMM::fitme(Yr[,i]~offset(l_ctot)+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2)
        print(bGLMMsA[[colnames(Yr)[i]]])
    }
    saveRDS(bGLMMsA, file=paste0(comb.out.dir, 'bGLMMsA.RDS'))

    #new code to extract dispersion estimates
    xab=rep(NA,length(bGLMMsA))
    for( i in 1:length(xab) ){
        test=bGLMMsA[[i]]
        #jfc, there must be a better way to extract this 
        a=capture.output(summary(test))[5]
        xab[i]=as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", a))
        print(xab[i])
    }
    names(xab)=names(bGLMMsA)
    saveRDS(xab, file=paste0(comb.out.dir, 'bGLMMsA_theta.RDS'))
    #system.time({f2=fitme(Yr[,3]~l_ctot+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='REML', data=data2, family=negbin2) })

    xab=readRDS(paste0(comb.out.dir, 'bGLMMsA_theta.RDS'))
    # code to process results of repeated measures mixed model with additional random additive effect  of genome
    bGLMMsA=readRDS(paste0(comb.out.dir, 'bGLMMsA.RDS'))

    ssH2a=matrix(0,length(bGLMMsA),6)
    #system.time({f3=fitme(Yr[,3]~offset(l_ctot)+expt+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2) })
    for( i in 1:length(bGLMMsA)) {    
        f3=bGLMMsA[[i]]
        lmbda=mean(exp(predict(f3)))
        #sigma for spaMM object isn't the dispersion estimate ... confusingly
        #sgma=sigma(f3)
        sgma= xab[names(bGLMMsA)[i]]
        #sgma=
        sA=VarCorr(f3)[1,3]
        sR=VarCorr(f3)[2,3]
        sC=VarCorr(f3)[3,3]
        ICC.A=sA/(sA+sC+sR+log(1+(1/lmbda)+(1/sgma)))
            #trigamma(1/((1/lmbda)+(1/sgma)))) #log(1+(1/lmbda)+(1/sgma)))
        ICC.R=sR/(sA+sC+sR+log(1+(1/lmbda)+(1/sgma)))
        ICC.C=sC/(sA+sC+sR+log(1+(1/lmbda)+(1/sgma)))
        ssH2a[i,1]=sA
        ssH2a[i,2]=sR
        ssH2a[i,3]=sC
        ssH2a[i,4]=ICC.A
        ssH2a[i,5]=ICC.R
        ssH2a[i,6]=ICC.C
    }
    rownames(ssH2a)=names(bGLMMsA)
    colnames(ssH2a)=c('AVar', 'RVar', 'ccVar', 'ICC.A', 'ICC.R', 'ICC.cc') #'ICC.cc', 'ICC.H2', 'ccVar', 'H2Var', 'CC.q', 'H2.q', 'tot.umis')

    saveRDS(ssH2a, file=paste0(comb.out.dir, 'scICC_3VC.RDS'))
    return(ssH2a)
}
