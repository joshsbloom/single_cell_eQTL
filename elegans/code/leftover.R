library(MASS)
library(pscl)
library(Rfast2)






  qm=qm2
    #wthl2=wthl[[tt]]
    #wthl2$chr=factor(wthl2$seqnames,levels=rev(c("X","V","IV","III","II","I")))
   


cPf$tpos=transcript.data$start_position[match(cPf$transcript, transcript.data$wbps_gene_id)]
cPf$CI.l=tstrsplit(cPf$CI.l, '_', type.convert=T)[[2]]
cPf$CI.r=tstrsplit(cPf$CI.r, '_', type.convert=T)[[2]]
cPf=cPf[cPf$tchr %in% c("X","V","IV","III","II","I"),]
cPf$tchr=factor(cPf$tchr,levels=c("X","V","IV","III","II","I"))
cPf$chr=factor(cPf$chr,levels=rev(c("X","V","IV","III","II","I")))

ggplot(cPf,aes(x=pos,y=tpos,alpha=-log10(FDR+1e-6)/6))+geom_point()+
    geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos,alpha=-log10(cPf$FDR+1e-6)/6)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none') +
    facet_grid(tchr~chr,scales="free")+theme_classic()


#tissueQTLPeaksNB
qm2=rbindlist(tissueQTLPeaksNB[-19],idcol='tissue')
tkk=names(which(sapply(tissueQTLPeaksNB[-19], function(x) sum(x$FDR<.05))>50))
qm2=qm2[qm2$tissue %in% tkk,]
#qm2=qm2[qm2$FDR<.05,]

qm2= tcPsig #rcp1f
qm2=qm2[  qm2$FDR< .05  & qm2$FDR.negbin<.05,] #qm2$FDR.negbin<.05 &
qm2$alpha=(-log10(qm2$FDR+1e-6)/6)
qm2$tpos=transcript.data$start_position[match(qm2$transcript, transcript.data$wbps_gene_id)]
qm2$CI.l=tstrsplit(qm2$CI.l, '_', type.convert=T)[[2]]
qm2$CI.r=tstrsplit(qm2$CI.r, '_', type.convert=T)[[2]]
qm2=qm2[qm2$tchr %in% c("X","V","IV","III","II","I"),]
qm2$tchr=factor(qm2$tchr,levels=c("X","V","IV","III","II","I"))
qm2$chr=factor(qm2$chr,levels=rev(c("X","V","IV","III","II","I")))
qm2=split(qm2, qm2$tissue)

wthl=lapply(withinTissueHotspots, function(x) {
           data.frame(markerGR[match(x[[1]], colnames(Gr)),])
  })
#wthl=rbindlist(lapply(wthl, data.frame),idcol='tissue')

test=lapply(af.tissue[tkk], function(x) {
           y=tstrsplit(names(x), '_', type.convert=T);
           data.frame(chr=y[[1]],pos=y[[2]], af=x) })
af.comb=rbindlist(test, idcol='tissue')
ggplot(af.comb, aes(x=pos, y=(1-af)-.5, color=tissue))+geom_line(size=1.25)+facet_grid(~chr)+ylab('allele freq. (obs-exp)')+
    geom_hline(yintercept=0,linetype="dotted")+
    theme_classic() 
ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/','AFdistortion.png'),width=15,height=7)

tot.count.LOD
tcL.comb=rbindlist(lapply(tot.count.LOD[names(which(sapply(tot.count.LOD, max)>4))], function(x) {
             y=tstrsplit(colnames(x), '_', type.convert=T);
           data.frame(chr=y[[1]],pos=y[[2]], LOD=x[1,]) }),idcol='tissue')
ggplot(tcL.comb, aes(x=pos, y=LOD, color=tissue))+geom_line(size=1.25)+facet_grid(~chr)+ylab('LOD for log2(total counts per tissue)')+
    theme_classic() 
ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/','LOD_tcount.png'),width=15,height=7)
ggplot(tcL.comb, aes(x=pos, y=LOD, color=tissue))+geom_line(size=1.25)+facet_grid(~chr)+ylab('LOD for log2(total counts per tissue)')+
    theme_classic()+ylim(0,20)
ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/','LOD_tcount_zoom.png'),width=15,height=7)



qm3=data.frame(qm2)
ggplot(qm3,aes(x=peakGpos,y=tGpos, alpha=alpha))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+  #facet_wrap(~tissue)+
    #geom_vline(aes(xintercept=ghot),data=wthl,color='red', size=1.25, alpha=.4)+
    geom_hline(yintercept=gcoord.key,color='darkblue', size=1)+
    geom_vline(xintercept=gcoord.key,color='darkblue', size=1)+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

library(data.table)
tkk=names(which(sapply(tissueQTLPeaksNB, function(x) sum(x$FDR<.05))>50))
wthsig=lapply(withinTissueHotspots, function(x) rbindlist(x[[2]]))
wthsig=wthsig[which(sapply(wthsig, function(x) sum(x$peaks))>0)]
wthsig$joint=rbindlist(jointHotspots[[2]])
hmerged=rbindlist(wthsig, idcol='tissue')

hmerged$tissue=relevel(as.factor(hmerged$tissue), ref="joint")
ggplot(hmerged, aes(x=values, y=lengths, color=peaks))+geom_col()+
    facet_grid(tissue~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('5cm bins across genome')


    ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/', 'hotspots.png'),width=12,height=10)







    scale_x_continuous(limits=c(0, max(gcoord.key)))+
    geom_hline(yintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=markerGR$gcoord[match(jointHotspots[[1]], colnames(Gr))], color='red', alpha=.9)+
    theme_classic()+
      theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

ggplot(cPf,aes(x=peakGpos,y=tGpos,alpha=alpha))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none') +
    scale_x_continuous(limits=c(0, max(gcoord.key)))+
    geom_hline(yintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=markerGR$gcoord[match(jointHotspots[[1]], colnames(Gr))], color='red', alpha=.9)+
    theme_classic()+
      theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())


jpsig=(jPs[jPs$FDR<.1,])
jpsigt=jpsig[jpsig$tchr!=jpsig$chr,]
hist(jpsigt$peakGpos,breaks=500)
abline(v=markerGR$gcoord[match(jointHotspots[[1]], colnames(Gr))],col='red')




#tissueQTLPeaksNB
qm2=rbindlist(tissueQTLPeaksNB,idcol='tissue')

#qm$peakGpos=markerGR$gcoord[match(qm$peak.marker, colnames(G))]
#qm$tGpos=transcript.data$gcoord[match(qm$transcript, transcript.data$wormbase_gene)]
#qm$peakLGpos=markerGR$gcoord[match(qm$CI.l, colnames(G))]
#qm$peakRGpos=markerGR$gcoord[match(qm$CI.r, colnames(G))]

#qm$markerpos=match(qm$peak.marker, colnames(G))
#qm$tpos=match(qm$transcript, transcript.data$wormbase_gene)
tkk=names(which(sapply(tissueQTLPeaksNB, function(x) sum(x$FDR<.05))>50))
qm2=qm2[qm2$tissue %in% tkk,]
#qm2=qm2[qm2$FDR<.05,]
qm2=qm2[qm2$FDR.negbin<.05 & qm2$FDR<.05,]
qm2$alpha=(-log10(qm2$FDR+1e-6)/6)
#qm2$alpha=-log10(qm2$FDR.negbin+1e-6)/10

wthl=lapply(withinTissueHotspots, function(x) {
           markerGR$gcoord[match(x[[1]], colnames(Gr))]
  })
wthl=do.call('c', wthl)
wthl=data.frame(tissue=gsub('\\..*', '', names(wthl)), ghot=wthl)

ggplot(qm2,aes(x=peakGpos,y=tGpos, alpha=alpha))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+  facet_wrap(~tissue)+
    geom_vline(aes(xintercept=ghot),data=wthl,color='red', size=1.25, alpha=.4)+
    geom_hline(yintercept=gcoord.key,color='darkblue', size=1)+
    geom_vline(xintercept=gcoord.key,color='darkblue', size=1)+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())


length(unique(do.call('c', lapply(tissueCisBetas,function(x) x$transcript[x$FDR.negbin<.05]))))
sum(cPs$FDR<.05)
usig=colnames(Yr)
cagg.betas=matrix(NA, length(tissueCisBetas), length(usig))
rownames(cagg.betas)=names(tissueCisBetas)
colnames(cagg.betas)=usig
for(tt in names(tissueCisBetas)){
    cagg.betas[tt,tissueCisBetas[[tt]]$transcript]=tissueCisBetas[[tt]]$negbin.beta
}

corGcis=cor(cagg.betas, use='pairwise.complete.obs')
corTcis=cor(t(cagg.betas), use='pairwise.complete.obs')


# left of 3 for Neuron
# first 2 on 5 for Intestine
# 


# 













cpss=cPs$transcript[cPs$FDR<.01]

disps=sapply(tissueNBNulls, function(x) sapply(x, function(y) y$info[4]))
disps= sapply(disps, function(x) {
                  y=x; 
                  names(y)=gsub('.dispersion', '',names(x)); 
                  return(y);})
udispn=unique(do.call('c', sapply(disps, names)))
disp.mat=matrix(NA,length(udispn),19)
rownames(disp.mat)=udispn
colnames(disp.mat)=names(disps)
for(tt in names(disps)){
    disp.mat[names(disps[[tt]]), tt]=disps[[tt]]
 }


cPf=tcP[tcP$FDR<.05,]
cPf$alpha=(-log10(cPf$FDR+1e-6)/6)


cPf=trcp1f[trcp1f$FDR.negbin<.05,]














    #scresid1=svd(scale(nbR))
    #plot(scresid1$d^2/sum(scresid1$d^2), xlim=c(0,20))
   
    #PCs=scale((scresid1$u[,1:1000]))
    #nt=as.integer(nrow(PCs))
    #rff=crossprod(PCs,scale(Gsub))/(nt-1)
    #LODr=-nt*log(1-rff^2)/(2*log(10))

    Yfe=nbR.R       # scale(nbR) #Rfast::colRanks(nbR, parallel=T)) #)#Yresid2.list$nbR)
    colnames(Yfe)=colnames(nbR)
    rownames(Yfe)=rownames(nbR)
    Gf=scale(Gsub) #Rfast::colRanks(Gsub, parallel=T))
    colnames(Gf)=colnames(Gsub)
    rownames(Gf)=rownames(Gsub)
   
    # Gsubr=(residuals(lm(lGreduced[tk,]~log(covs$total)+covs$Batch)))
   # GsubrE=exp(Gsubr)/(1+exp(Gsubr))
    #Gsubr=(residuals(lm(Gsub~log(covs$total)+covs$Batch)))
    #Gf=apply(Gsubr, 2, rankNorm)
    #Gsubr.var=colVars(Gsubr)
   # Gf=scale(Gsubr)

    #dataC=data.frame(cbind(mmp1, Gf))
     #library(mpath)
    #gn='WBGene00020612'
    #test= glmreg(Yk[,gn]~., data=dataC, family="negbin", theta=cisNB[[gn]]$theta)

    nt=as.integer(nrow(Yfe))
    rff=crossprod(Yfe,Gf)/(nt-1)
    
    # don't retain standardized betas within tissue given difficulty 
    # of interpreting across tissues 
    #saveRDS(rff, paste0(dout, 'tissue_std_beta.RDS'))
    LODr=-nt*log(1-rff^2)/(2*log(10))
    
    ##saveRDS(fl(LODr), paste0(dout, 'tissue_LOD.RDS'))
    # and construct permutation nulls 
    nperm=50
    LODperm=replicate(nperm, {
                    nind=sample(1:nt)
                    pme=crossprod(Yfe[nind,],Gf)/(nt-1)
                    return((-nt*log(1-pme^2)/(2*log(10)))) 
                    })
    #lpmax=apply(LODperm, c(1,3), max)
    
   
    tchr=transcript.data$chromosome_name[match(rownames(LODr),transcript.data$wbps_gene_id)]
    mchr=as.character(seqnames(mGRs))
    eg=expand.grid(tchr, mchr)
    be=which(eg$Var1==eg$Var2)

    #for(column in 1:ncol(emat)){emat[rownames(emat)==colnames(emat)[i],column]=1; print(column); }
    LODr[be]=0
    for(i in 1:nperm) {     LODperm[,,i][be]=0    }
    
    # if retaining permutations 
    LODmatrix[[uc]]=LODr
    LPmatrix[[uc]]=LODperm


    #m1=which(mchr=='I')
    #pmean=apply(LODperm[,m1,],c(2,3), mean)
    #plot(colMeans(LODr[,m1]))

    mGRs=markerGR[match(colnames(Gf), colnames(Gr))]
    twgFDR= getFDRfx(LODr,LODperm) #, cisMarkers) 
   
       
   
    tcP=getGlobalPeaks(rff, LODr,twgFDR, mGRs,  transcript.data, chroms) 
    tcPsig=fitNBTransModel(tcP, Yk,Gsub, cisNB, mmp1, cMsubset, .5) 
      
    t.hotspots = getHotspots(tcPsig, gmapd, fdr.thresh=.1, pFDR=T)
    ggplot(data.frame(rbindlist(t.hotspots[[2]])), aes(x=values, y=lengths, color=peaks))+geom_col()+
        facet_grid(~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('10cm bins across genome')
    x11()
    t.hotspotsF = getHotspots(tcPsig, gmapd, fdr.thresh=.1, pFDR=F)
    ggplot(data.frame(rbindlist(t.hotspotsF[[2]])), aes(x=values, y=lengths, color=peaks))+geom_col()+
        facet_grid(~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('10cm bins across genome')

    # where is the na coming from?
    qm=na.omit(tcPsig[ tcPsig$FDR.negbin<(10^-3),])
    ggplot(qm,aes(x=pos,y=tpos,alpha=-log10(FDR.negbin+1e-10)))+geom_point()+
        geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) + #, color ="grey70", alpha=.1))+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') +  
        #geom_vline(aes(xintercept=start),data=wthl2,color='red',alpha=.75)+
        facet_grid(tchr~chr,scales="free")+theme_classic()  



sG=svd(Gf)
sum(sqrt(sG$d))^2/sum(sG$d)
sN=svd(nbR)
sum(sqrt(sN$d))^2/sum(sN$d)







    #af.tissue[[uc]]=colSums(Gsub)/nrow(Gsub)
    cprod=crossprod(scale(log2(covs$total)),scale(Gsub))/(nrow(Yk)-1)
    tot.count.LOD[[uc]]=-nrow(Yk)*log(1-cprod^2)/(2*log(10))
    #======================================================================================

  
   
    #"V_536997"
    #test=betareg(Gsub[,174]~log(covs$total)+covs$Batch)
    GsubA=((lm(lGreduced[tk,]~log(covs$total)+covs$Batch)))
    
    #prevent numerical issues    
    Gsub[Gsub==1]=.9999999
    Gsub[Gsub==0]=.0000001 

    Gsubr=(residuals(lm(lGreduced[tk,]~log(covs$total)+covs$Batch)))
    GsubrE=exp(Gsubr)/(1+exp(Gsubr))
    Gsubr.var=colVars(Gsubr)
    Gsubr=standardise2(GsubrE)
    #======================================================================================
   
    #run once =============================================================================
    # goal here is to estimate dispersion within each tissue for each transcript
    #nbM0=negbinNull(colnames(Yk),Yk,mmp1)
    #saveRDS(nbM0, paste0(dout,'nb_null.RDS'))
    #tissueNBNulls[[uc]]=nbM0
    nbM0=tissueNBNulls[[uc]]
 
    ## recalculate betas at joint hotspots, requires defining jointHotspots[[1]] ==========
    ## which happens outside this loop
    #jnbHotspots=negbinRefitHotspots(jointHotspots[[1]], Yk, Gsub,mmp1,nbM0)
    #saveRDS(jnbHotspots,paste0(dout, 'tissue_joint_hotspotBetas.RDS'))
    #jointTissueHotspotBetas[[uc]]=jnbHotspots
    ##=====================================================================================

    # for OLS Log2(Y+1) model get standardized genotypes ==================================    
    #  Gf=standardise2(Gsub)
    # colnames(Gf)=colnames(Gsub)
    # rownames(Gf)=rownames(Yk)
    #======================================================================================

    #nbNullRefit=negbinNullDev(tnames, Yk, mmp1,nbM0, smoothed.disp)
    Yresid2.list=negbinNullDev(colnames(Yk), Yk, mmp1,nbM0, covs)
    #Yresid=negbinNullDev(tnames, Yk, mmp1,nbM0)
    #Yfe=standardise2(Yresid[rownames(Gsubr),])

    Yfe=standardise2(Yresid2.list$nbR)
    Gf=Gsubr
    nt=as.integer(nrow(Yfe))
    rff=crossprod(Yfe,Gf)/(nt-1)
    
    # don't retain standardized betas within tissue given difficulty 
    # of interpreting across tissues 
    #saveRDS(rff, paste0(dout, 'tissue_std_beta.RDS'))
    LODr=-nt*log(1-rff^2)/(2*log(10))
    LODmatrix[[uc]]=LODr
    
    ##saveRDS(fl(LODr), paste0(dout, 'tissue_LOD.RDS'))
    # and construct permutation nulls 
    nperm=10
    LODperm=replicate(nperm, {
                    nind=sample(1:nt)
                    pme=crossprod(Yfe[nind,],Gf)/(nt-1)
                    return((-nt*log(1-pme^2)/(2*log(10)))) 
                    })
    LPmatrix[[uc]]=LODperm
    ### faster to keep in memory than save and reload 
    ###saveRDS(apply(LODperm,3,fl), paste0(dout, 'tissue_LOD_perm.RDS'))
    #======================================================================================


    # return functions to map LOD to FDR for global and cis-only tests=====================
    twgFDR= getFDRfx(LODr,LODperm, cisMarkers) 
    qtl.fdrs[[uc]]=twgFDR
    #saveRDS(twgFDR, paste0(dout, 'FDRfx.RDS'))
    #twgFDR=qtl.fdrs[[uc]]
    #======================================================================================
    
    # test the closest marker to each transcript ===========================================
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(rff)),]
    cisBetasTis=rff[cbind(cMsubset$transcript, cMsubset$marker)]
    names(cisBetasTis)=cMsubset$transcript
    max.obsLODc=LODr[cbind(cMsubset$transcript, cMsubset$marker)]

    cisBetasTis=data.frame(transcript=names(cisBetasTis), 
                           cisMarker=cMsubset$marker,
                           sbeta=cisBetasTis,
                           LOD=max.obsLODc,
                           FDR=twgFDR[[4]](max.obsLODc), stringsAsFactors=F)
    cisBetasTis$FDR[is.na(cisBetasTis$FDR)]=1
    cisBetasTis$FDR[cisBetasTis$FDR>1]=1
    
    #should refit cis betas here!!
    cisNB=cisNegbinFit(cMsubset,Yk,Gsub,mmp1,nbM0,Yresid2.list)
    cisNB=rbindlist(cisNB,fill=T,idcol='transcript')
    cisBetasTis=merge(cisBetasTis, cisNB, by='transcript', all.x=T,sort=F)
    cisBetasTis$FDR.negbin=p.adjust(cisBetasTis$negbin.p, method='fdr')
    tissueCisBetas[[uc]]=cisBetasTis
    #saveRDS(cisBetasTis, paste0(dout, 'tissue_cisBetas.RDS'))
    #========================================================================================

    # get peaks from global analysis within each tissue ===================================
    mGRs=markerGR[match(colnames(Gf), colnames(Gr))]
    tcP=getGlobalPeaks(rff, LODr,twgFDR, mGRs,
                        transcript.data, chroms) 

    tissueQTLPeaks[[uc]]=tcP
    # saveRDS(tcP, paste0(dout,'Peaks.RDS')) 
    #tcP=tissueQTLPeaks[[uc]]  
    #======================================================================================


    # fit negative binomial regression model for each global peak if FDR < .2 =============
    tcP_nb=negbinRefit(tcP, Yk,Gsub, mmp1, nbM0, fdr.thresh=.1,Yresid2.list)
    trcp1=rbindlist(tcP_nb,fill=T)
    trcp1f=trcp1[!is.na(trcp1$negbin.p),]
    trcp1f$FDR.negbin=p.adjust(trcp1f$negbin.p, method='fdr')
    # saveRDS(trcp1f, paste0(dout,'PeaksNB.RDS')) 
    tissueQTLPeaksNB[[uc]]=trcp1f #tcP
    #trcp1f= tissueQTLPeaksNB[[uc]] 
    #======================================================================================
    
    # analysis of hotspots within each tissue =============================================
    t.hotspots = getHotspots(trcp1f, gmapd, fdr.thresh=.05, pFDR=F)
    #t.hotspots= withinTissueHotspots[[uc]]
    withinTissueHotspots[[uc]]=t.hotspots
    # saveRDS(t.hotspots,paste0(dout, 'within_tissue_hotspots.RDS'))
    
    # calc negative binomial model and within-tissue hotspots 
    if( length(t.hotspots[[1]])>0 ) {
      #uncomment  wtHotspots=negbinRefitHotspots(t.hotspots[[1]], Yk, Gsub,mmp1,nbM0)
        #withinTissueHB=rff[,tissue.hotspots[[1]]]
        #saveRDS(wtHotspots,paste0(dout, 'within_tissue_hotspot_betas.RDS'))
     #uncomment   withinTissueHotspotBetas[[uc]]=wtHotspots
    }
    #=======================================================================================


   }


      # OLS Log2(Y+1) model, get standardized residuals correcting for total reads ==========
    # and batch effects
    # presumably run once only
    BOLS=lm(log2(Yk+1)~log2(covs$total)+covs$Batch)
    Yresid=residuals(BOLS)
    rownames(Yresid)=rownames(Yk)
    colnames(Yresid)=colnames(Yk)
    rm(BOLS)
    #residual.matrices[[uc]]=Yresid
    gc()
    #so slow ...
    Yfe=standardise2(Yresid)
    rownames(Yfe)=rownames(Yresid)
    colnames(Yfe)=colnames(Yresid)
    # transcripts that are singular after the regression 
    # (no variance = nothing to map on)
    # but note that next two lines are legacy and aren't necessary 
    shitT=unique(which(is.na(Yfe), arr.ind=T)[,2])
    if(length(shitT)>0) {    Yfe=Yfe[,-shitT] } 
    #======================================================================================


   
    # Fast matrix-based mapping ===========================================================
        # then map within each classification

    # contrast with deviance residuals
        # remove dud transcripts here (To DO )))

    #rownames(Yfe)=rownames(Yresid)
    #colnames(Yfe)=colnames(Yresid)
    
    #diagnostics
    #cor_w_total=cor(log2(Yk+1),log2(covs$total))
    #plot(log2(colSums(Yk)), cor_w_total, xlab='log2 total expression for transcript in BWM', ylab='r (log(transcript) v log(total)', col="#00000055")
    tgns=colnames(Yk)[which(cor_w_total>.6)]
    tgn=tgns[3]

    tgn="WBGene00017210" #tgns[30]
    par(mfrow=c(2,1))
    plot(LODr[tgn,], ylab='LOD', main='OLS LOD', ylim=c(0,60))

    nmodelf=nbM0[[tgn]]
    gll=rep(NA,ncol(Gsub))
    names(gll)=colnames(Gsub)
    for(p in colnames(Gsub)[seq(1,ncol(Gsub),10)]){
        print(p)
       gll[p]=logLik(fastglm(cbind(mmp1,Gsub[,p]), Yk[,tgn], family=negative.binomial(theta=nmodelf$info[4],link='log')))
    }
        plot(-2*(nmodelf$info[3]-gll)/2.3, main='glm LOD', ylab='LOD') 
      
    tsplit=split(trcp1f, trcp1f$transcript)
    tsplitc=sapply(tsplit,nrow)



    
       
     tcP.OLS=trcp1f
     H.OLS=getHotspots(tcP.OLS, gmapd, fdr.thresh=.05, pFDR=F)
     dff=rbindlist(H.OLS[[2]],idcol='chrom')
    ggplot(dff, aes(x=values, y=lengths, color=peaks))+geom_col()+
     facet_grid(~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('5cm bins across genome')


     H.glm=getHotspots(tcP, gmapd, fdr.thresh=.05, pFDR=T)
     dff=rbindlist(H.glm[[2]],idcol='chrom')
    ggplot(dff, aes(x=values, y=lengths, color=peaks))+geom_col()+
     facet_grid(~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('5cm bins across genome')













                     apply(x,1,min, na.rm=T) })
    minpO=apply(p, 1, min, na.rm=T)
    prange=range(minp0)


    rvec=-log10(apply(p,1,min))
    nvec=(-log10(apply(pp,1, min)))

    Meff=37
    plot(which(p<.001, arr.ind=T))
    #gchr=tstrsplit(colnames(Gsub), '_', type.convert=T)[[1]]
    #minpc=apply(p, 1, function(x) { sapply(split(x, gchr), min) } )
    #mps=apply(minpc, 2, function(h) {   -cumsum(log(1-sort(h*Meff)))/seq_along(h)     })
    #sum(sort(sapply(mps, function(x) sum(x<.05, na.rm=T))))


    stopCluster(cl)






    return(testme)
}
          


    LRSs=foreach(gn=names(cisNB)[1:128], .combine=rbind,
                 .export=c("cisNB", "Yk", "Gsub", "mmp1", "nbLL" ),
                 .packages=c('fastglm')) %dopar% {
       print(match(gn, names(cisNB)))
       theta.est=cisNB[[gn]]$theta
       fnbrN=fastglmPure(mmp1, Yk[,gn],
                         family=negative.binomial(theta=theta.est,link='log'), maxit=100)
       nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)

       LRS=rep(NA,ncol(Gsub))
       names(LRS)=colnames(Gsub)     
       for(gx in colnames(Gsub)){
           XX=cbind(mmp1, Gsub[,gx]) 
           fnbrF=fastglmPure(XX, Yk[,gn],family=negative.binomial(theta=theta.est,link='log'), maxit=100)
           LRS[gx]=-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))
        }
       return(LRS)
    }

    rownames(LRSs)=names(cisNB)

    #LRSm=do.call('rbind', LRSs)
    rownames(LRSm)=names(cisNB)
    return(LRSm)
}


#nmLLik=cisNB[[gn]]$LLik
#class(fnbrN)='glm'
#nmLLik=as.vector(logLik(fnbrN0))


    cGsub=cor(Gsub)    
    Meff=getMeff_Li_and_Ji(cGsub)

    p=pchisq( LRSm,1, lower.tail=F)
    gchr=tstrsplit(colnames(Gsub), '_', type.convert=T)[[1]]
    minpc=apply(p, 1, function(x) { sapply(split(x, gchr), min) } )
    mps=apply(minpc, 2, function(h) {   -cumsum(log(1-sort(h*Meff)))/seq_along(h)     })
    sum(sort(sapply(mps, function(x) sum(x<.05, na.rm=T))))
    
    nlp=-log10(p)

    apvec=sort(as.vector(p), decreasing=F)

    qtransf=getFDRmeff(apvec,q=.05,Meff)


    
    mpc=svd(Gsubr)
    top10eig=mpc$u[,1:5]
    mmp2=scale(cbind(mmp1, top10eig))
    nbM2=negbinNullP(colnames(Yk), Yk, scale(mmp1))
    
    #transcripts with NA dispersions tend to be very low expressed 
    str(Yk[,which(disps2>100)])
    
    which(disps2>100)

    #nbM1=negbinNull(colnames(Yk),Yk,cbind(mmp1,top10eig))
    disps2=sapply(nbM1, function(x)x$info[4])
    names(disps2)=gsub('\\.dispersion', '', names(disps2))

        

    nbM1=negbinNullP(colnames(Yk)[1:100],Yk,mmp2)

    nbM2=negbinNullP(colnames(Yk)[1:100],Yk,cbind(mmp1, top10eig[,10]))

    # logic from this ###################################################################
    #https://github.com/ChristophH/sctransform
    disps=sapply(nbM0, function(x)x$info[4])
    names(disps)=gsub('\\.dispersion', '', names(disps))
#
    bad.fits=which(log10(disps)>2)
    nbbest=0
    iter=0
    for(gn in colnames(Yk)){
        iter=iter+1
        print(gn)
        test=mgcv::gam(Yk[,gn]~cbind(mmp1[,-1]), family=mgcv::nb , method="ML")
        test$family$getTheta(T)
        test3=glm.nb(Yk[,gn]~cbind(mmp1[,-1],top10eig))
        test$family$getTheta(T)


        test2=fastglm(cbind(mmp1),Yk[,gn],family=poisson(link='log'), method=3)
        pt2=predict(test2,cbind(mmp1))
        theta.ml(Yk[,gn], exp(pt2))

        if(BIC(test)<BIC(test2)) {nbbest=nbbest+1 }
        print(paste(BIC(test), BIC(test2), BIC(test)<BIC(test2)))
        print(nbbest/iter)
    }
#        test=glm.nb(Yk[,gn]~mmp1[,-1]) #, family=mgcv::nb , method="ML")
#        print(test$theta)
#        #print(test$family$getTheta(T))
#    }
#
    mexp=colMeans(Yk)
    l10disps=log10(disps2)
    l10mexp=log10(mexp)
    l10dispsf=l10disps[l10disps>-4 & l10disps<4]
    l10mexpf=l10mexp[l10disps>-4 & l10disps<4]
    bw=bw.SJ(na.omit(l10mexpf))*3
    ksy=ksmooth(x=l10mexpf, y=l10dispsf, x.points=l10mexp, bandwidth=bw, kernel='normal')
    m.to.theta2=approxfun(10^ksy$x, 10^ksy$y, ties=mean, rule=2)
    #=###################################################################################



#investigation of genotype informative snps per effect on mapping 

    for(cc in chroms) {
        print(cc)
        Ns=N2.counts[grep(paste0('^', cc, '_'), rownames(N2.counts)),covs$barcode]
        Cs=CB.counts[grep(paste0('^', cc, '_'), rownames(N2.counts)),covs$barcode]
        inf3=sum(colSums(Cs)>3 & colSums(Ns)>3)
        inf0N=sum(colSums(Ns)==0)
        inf0C=sum(colSums(Cs)==0)
        print(paste('tot cells=', ncol(Ns)))
        print(paste('more than 3 informative reads, both genotypes',inf3))
        print(paste('cells with no informative N2 reads', inf0N))
        print(paste('cells with no informative CB reads', inf0C))
        #Ns=N2.counts[,covs$barcode]
        #Cs=CB.counts[,covs$barcode]
        png(file=paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/intestine', cc, '.png'), width=1024, height=1024)
        plot(jitter(log2(colSums(Ns)),40), main=cc,
         jitter(log2(colSums(Cs)),40), xlim=c(0,10), ylim=c(0,10), col="#00000044", 
         xlab='jittered log2 informative N2 sites per cell - Intestine',
         ylab='jittered log2 informative CB sites per cell - Intestine')
        dev.off()
        
    }
   Ns=N2.counts[,covs$barcode]
   Cs=CB.counts[, covs$barcode]
   plot(jitter(log2(colSums(Ns)),40),         
        jitter(log2(colSums(Cs)),40), col="#00000044", xlab='log2(N informative sites)', ylab='log2(C informative sites)', main='intestine')
   abline(h=4.17)
   abline(v=4.17)

   tkf1=which((colSums(Cs)>18 & colSums(Ns)>18))
   tk=tk[tkf1]
   covs=covs[tkf1,]



Ns=N2.counts[grep('^X_', rownames(N2.counts)),covs$barcode]
Cs=CB.counts[grep('^X_', rownames(N2.counts)),covs$barcode]

ggplot(data.frame(Ns=log2(colSums(Ns)+1), Cs=log2(colSums(Cs)+1)), aes(Ns,Cs))+geom_bin2d()

plot(jitter(log2(colSums(Ns)),40), 
     jitter(log2(colSums(Cs)),40), xlim=c(0,10), ylim=c(0,10), col="#00000044", 
     xlab='jittered log2 informative N2 sites per cell - Hypodermis',
     ylab='jittered log2 informative CB sites per cell - Hypodermis',
     sub='779/13219 cells have >3 N2 geno informative reads and >3 CB geno informative reads'

)


suspect.set=(trcp1f[trcp1f$chrom=='I' &trcp1f$FDR<.01 &trcp1f$tchr!='I',])





Gsub









mmp3=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch))

mmp0=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch-1))
#mmp1=cbind(mmp0, Gf[,hotspots$hotspot.markers])
    
mmp00=model.matrix(lm(log2(Yk[,1]+1)~covs$Batch-1))

     # first use lasso to fit
    mmp2=cbind(mmp1,Gf[,ptrq[[gn]]$peak.marker])
    lasso <- cv.glmnet(
            y=Yk[,gn],
            x=mmp2,
            family="poisson",
            alpha=0 # <- lasso
        )
        pred <- as.vector(exp(predict(lasso, newx=mmp2))) #, newoffset=offset)))
        theta <- theta.ml(y=Yk[,gn], mu=pred)

    for(gn in names(ptrq)) {
        print(gn)
        print(ptrq[[gn]])

        #mmp1=cbind(mmp00, Gsub[,ptrq[[gn]]$peak.marker])
        #testz=(zeroinfl(Yk[,gn]~offset(log(covs$total))+covs$Batch+Gsub[,ptrq[[gn]]$peak.marker], dist='negbin'))
        testnb=glm.nb(Yk[,gn]~offset(log(covs$total))+covs$Batch+Gsub[,as.character(ptrq[[gn]]$peak.marker)]) 
        testnb2=glm.nb(Yk[,gn]~(log(covs$total))+covs$Batch+Gsub[,as.character(ptrq[[gn]]$peak.marker)] ) 

        print(summary(testnb))
        testlm=lm(log2(Yk[,gn]+1)~(log2(covs$total))+covs$Batch+Gsub[,as.character(ptrq[[gn]]$peak.marker)] -1)
        print(summary(testlm))
        gn="WBGene00000475"
    }

        testp=(glm(Yk[,gn]~offset(log(covs$total))+covs$Batch+Gsub[,ptrq[[gn]]$peak.marker],family='poisson'))
        mmp4=cbind(mmp3,scale(Gf[,as.character(ptrq[[gn]]$peak.marker)]))
        testnb2=(glm.nb(Yk[,gn]~mmp4))
        
        testnrn= negbin.reg(Yk[,gn], mmp3) #4[,-1])
        ds=seq(1,ncol(Gf),10)
        
        
        lll=rep(NA,ncol(Gf))
        llf=rep(NA,ncol(Gf))

        for(i in ds){
            print(i)
            mmp4=cbind(mmp3,Gf[,i])
            lll[i]=as.numeric(negbin.reg(Yk[,gn], mmp4[,-1], maxiters=500)$info[3])
        }
          plot(-2*(testnrn$info[3]-lll)/(2*log(10)), ylim=c(0,4))   

        print(ptrq[[gn]])
        print(summary(testnb))
        print(summary(lm(log2(Yk[,gn]+1)~log2(covs$total)+covs$Batch+Gsub[,ptrq[[gn]]$peak.marker])))
    }
        #test=(glm.nb(Yk[,gn]~(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker]))

        test=glm.nb(Yk[,gn]~offset(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker],family= negative.binomial(theta=.3969))

        test=(glmmTMB(Yk[,gn]~offset(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker], family='nbinom2'))

        test1=glm(Yk[,gn]~offset(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker],family=negative.binomial(theta=.3969))
        test1=glm(Yk[,gn]~(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker],family=negative.binomial(theta=.34))

        test1=glm(Yk[,gn]~offset(log(covs$total))+mmp00+Gsub[,ptrq[[gn]]$peak.marker],family='poisson')

    }

    test=glm.nb(Yk[,'WBGene00269355']~mmp) #, family='poisson')z


cisNB=list()
for(gn in cisMarkers$transcript){
    #gn=colnames(Yr)[1]
    print(gn)
    cMarker=cisMarkers$marker[cisMarkers$transcript==gn]
    g00=glm(Yr[,gn]~mm-1, family='poisson')
    tech.covar.n=length(coef(g00))

    theta.est=theta.ml(Yr[,gn], exp(predict(g00)), limit=100)
    
    ctest=tryCatch({
      g0=update(g00, ~., family=negative.binomial(theta=theta.est))
      g1=update(g0, ~.+Gr[,cMarker])
      cisBeta=coef(g1)[tech.covar.n+1] 
      cisLRS=-2*(logLik(g0)-logLik(g1))
      cisP=pchisq(cisLRS , 1, lower.tail=F)
      return(list(cisBeta=cisBeta, cisLRS=cisLRS, cisP=cisP,g0=g0))},
      error=function(e) {
           g1=update(g00, ~.+Gr[,cMarker])
           cisBeta=as.numeric(coef(g1)[tech.covar.n+1] )
           cisLRS=NA
           cisP=NA #pchisq(cisLRS , 1, lower.tail=F)
          return(list(cisBeta=as.numeric(cisBeta), cisLRS=cisLRS, cisP=cisP))

      })
   cisNB[[gn]]=c('Beta'=ctest$cisBeta, 'LRS'=ctest$cisLRS, 'P'=ctest$cisP)
   print(cisNB[[gn]])
   g0=ctest$g0
   print(gn %in%  names(cP_nb))
   if(is.null(g0)) {next;}

   if(gn %in%  names(cP_nb) ) {

        pmarker.vec=as.character(cP_nb[[gn]]$peak.marker)
        cP_nb[[gn]]$negbin.beta=rep(NA, length(pmarker.vec))
        cP_nb[[gn]]$negbin.LRS=rep(NA, length(pmarker.vec))
        cP_nb[[gn]]$negbin.p=rep(NA, length(pmarker.vec))
        for(pp in 1:length(pmarker.vec)) {
            p=pmarker.vec[pp]
            sout=tryCatch({
            g2=update(g0,~.+Gr[,p])
            return(list(  
                      Beta=coef(g2)[tech.covar.n+1],
                      LRS=-2*(logLik(g0)-logLik(g2)),
                      P=pchisq(LRS , 1, lower.tail=F)))
            },error=function(e){return(list(Beta=NA,LRS=NA,P=NA)) } )
            print(sout)
            cP_nb[[gn]]$negbin.beta[pp]=sout$Beta
            cP_nb[[gn]]$negbin.LRS[pp]=sout$LRS
            cP_nb[[gn]]$negbin.p[pp]=sout$P
        }
        print(cP_nb[[gn]])
    }

   print(match(gn,   cisMarkers$transcript))
}

cisNB=list()
for(gn in cisMarkers$transcript){
    #gn=colnames(Yr)[1]
   
    cMarker=cisMarkers$marker[cisMarkers$transcript==gn]
    #g00=glm(Yr[,gn]~mm-1, family='poisson')
    #theta.est=theta.ml(Yr[,gn], exp(predict(g0)), limit=50)
    #g0=update(g00, ~., family=negative.binomial(theta=theta.est))
    g0=negbin.reg(Yr[,gn],mm[,-1])
    
    tech.covar.n=length(g0$be)
  
    g1=negbin.reg(Yr[,gn],cbind(mm[,-1],Gr[,cMarker] ))
    g1=glm.nb((Yr[,gn]~cbind(mm[,-1],Gr[,cMarker] ))
    cisBeta=g1$be[tech.covar.n+1,1] 
    cisLRS=as.numeric(-2*(g0$info[3]-g1$info[3]))
    cisP=pchisq(cisLRS , 1, lower.tail=F)
    cisNB[[gn]]=c(cisBeta,cisLRS,cisP)
    print(cisNB[[gn]])
    if(gn %in%  names(cP_nb) ) {
        pmarker.vec=as.character(cP_nb[[gn]]$peak.marker)
        cP_nb[[gn]]$negbin.beta=rep(NA, length(pmarker.vec))
        cP_nb[[gn]]$negbin.LRS=rep(NA, length(pmarker.vec))
        cP_nb[[gn]]$negbin.p=rep(NA, length(pmarker.vec))
        for(pp in 1:length(pmarker.vec)) {
            p=pmarker.vec[pp]
            g2=negbin.reg(Yr[,gn],cbind(mm[,-1],Gr[,p] ))
            Beta=g2$be[tech.covar.n+1,1] 
            LRS=as.numeric(-2*(g0$info[3]-g2$info[3]))
            P=pchisq(LRS , 1, lower.tail=F)
            cP_nb[[gn]]$negbin.beta[pp]=Beta
            cP_nb[[gn]]$negbin.LRS[pp]=LRS
            cP_nb[[gn]]$negbin.p[pp]=P
        }
        print(cP_nb[[gn]])
    }

   
    
    
    print(match(gn,   cisMarkers$transcript))
}






        
        cP_1[[gn]]$negbin.beta=plb
        cP_1[[gn]]$negbin.LRS=-2*(nmodel-plr)
        cP_1[[gn]]$negbin.p=pchisq( cP_1[[gn]]$negbin.LRS,1, lower.tail=F)  #-2*(nmodel-plr)
    }
   
    

    #nmodel= negbin.reg(Yr[,gn], mm[,-1], maxiters=100)
    system.time({g1=glm(Yr[,gn]~mm+Gr[,cMarker]-1, family=negative.binomial(theta=theta.est))} ) # #'poisson'))
    #system.time({ g2=negbin.reg(Yr[,gn], cbind(mm[,-1],Gr[,cMarker]), maxiters=100)})

    
    #$info[3] #[,-1])    

}


#cP=readRDS('/data/single_cell_eQTL/elegans/results/20191217/joint/jointPeaks.RDS')
#refit
#cisBetasCount=rep(NA, length(cisMarkers$transcript))
#names(cisBetasCount)=cisMarkers$transcript
#tech.covar.n=ncol(mm)
#for(gn in cisMarkers$transcript){
#    print(gn)
    #g00=speedglm(Yr[,gn]~mm-1, , family=poisson(link='log'), fitted=T) #family='poisson')
    #theta.est=theta.ml(Yr[,gn], exp(predict(g00)), limit=100)
    #g01=speedglm(Yr[,gn]~mm-1, , family=negative.binomial(theta=theta.est))
    #tech.covar.n=length(coef(g00))
    #theta.est=theta.ml(Yr[,gn], exp(predict(g00)), limit=100)
#    cMarker=cisMarkers$marker[cisMarkers$transcript==gn]
#    cisBetasCount[gn]=speedglm(Yr[,gn]~ cbind(mm[,-1], Gr[,cMarker]), family=poisson(link='log'))$coef[tech.covar.n+1]
#}

#8452


#mm2=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$fine_tissue) #+joint.covariates$fine_tissue)
#BOLS=lm.fit(mm2, log2(Yr+1), intercept=F)
#BC=(BOLS$coefficients[-c(1:6),])
#rownames(BC)=gsub("joint.covariates.fine_tissue", "", rownames(BC))

#par(oma=c(10,10,10,10))
#heatmap.2(BC, col=redgreen(75), trace='none',scale='column')
#test=umap(t(BC))
#
#saveRDS(BC, file='/data/single_cell_eQTL/elegans/results/perFineTissueBetas.RDS')
#saveRDS(test, file='/data/single_cell_eQTL/elegans/results/perFineTissueUMAP.RDS')
#
#indexp=(apply(BC,1, function(x) (abs(x>.25))))
#tbetathresh=abs(BC)>.2

#plot(test$layout[,1], test$layout[,2], xlim=c(-10,10),ylim=c(-10,10), col=1+(tbetacnt>0))
#for(i in 1:ncol(indexp))  {
#plot(test$layout[,1], test$layout[,2], xlim=c(-10,10),ylim=c(-10,10), col=indexp[,i]+1, main=colnames(indexp)[i], xlab='umap on tissue betas 1', ylab='umap on tissue betas 2')
#readline()
#}

#indices of cis eQTL in 
#rows are transcripts, columns are markers
#cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
#------------------------------------------------------------------------------------------------------------
#sinfo=seqinfo(markerGR)
#sinfo@seqlengths=sapply(sapply(sapply(split(markerGR, seqnames(markerGR)),function(x) start(x)), max), as.integer)
#Gfit=residuals(lm.fit(mm, G, intercept=F))
#Gfit=standardise2(Gfit)
#rownames(Gfit)=rownames(G)
#colnames(Gfit)=colnames(G)

# 1D scan speedup 
#rt2=covar(Yre,G) #rossprod(Yre,G)/(nrow(Yre)-1)
jcov2=joint.covariates
jcov2$broad_tissue=as.factor(jcov2$broad_tissue)

joint.resids=matrix(NA,nrow(Yr),length(expressed.transcripts))
rownames(joint.resids)=rownames(Yr)
colnames(joint.resids)=expressed.transcripts

#nbNull=list()
for(gn in expressed.transcripts){
    print(gn)
    lookup=joint.covariates$broad_tissue %in% names(which(is.expressed20[,gn]))
    jck=jcov2[lookup,]
    jck=droplevels(jck)

    Ys=Yr[lookup,gn]
    if(length(levels(jck$broad_tissue))==1) {
        m0l=model.matrix(Ys~log2(jck$total)+jck$Batch)
    } else {
        m0l=model.matrix(Ys~log2(jck$total)+jck$Batch+jck$broad_tissue)
    }
    joint.resids[lookup,gn]=residuals(lm.fit(m0l,Ys))
    #nbNull[[gn]]=negbin.reg(Ys,m0l[,-1])
}
    
#test= fastglm(m0l, Ys, family=negative.binomial(theta=nm$info[4],link='log'))
#rvec=cor(residuals(test), Gr[lookup,])
#test2= fastglm(m0l, log2(Ys+1), family=gaussian())
#rvec2=cor(residuals(test2), Gr[lookup,])
#rvec=cor(residuals(test), Gr[lookup,])
#test= glm(Ys~m0l, family=negative.binomial(theta=nm$info[4],link='log'))
# t1= fastglm(m0l, log2(Ys+1), family=gaussian() ) #negative.binomial(theta=nm$info[4],link='log'))

ncount=apply(joint.resids,2,function(x) sum(!is.na(x)))

Yre=apply(joint.resids,2, function(x) {
                goodx=!is.na(x) #x!=0
                y=x[goodx]
                y=scale(y)
                x[goodx]=y
                return(x)
                })
Yre[is.na(Yre)]=0

Yre2=Yre
r2=crossprod(Yre2,G)/(ncount-1)
tt=(r2/sqrt(1-r2^2))*sqrt(ncount-2)
ptt=(2*pt(-abs(tt),df=ncount-2))














    #new poisson model ------------------------------------------
    mmp0=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch-1))
    mmp1=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch))

    #unstandardized genotypes 
    mmp=cbind(mmp0, Gsub[,hotspots$hotspot.markers])
       
    for(ty in 1:ncol(Yk)) {
       print(ty)
       glms[[uc]][[colnames(Yk)[ty]]]= tidy(glm(Yk[,ty]~mmp, family='poisson'))
       lms[[uc]][[colnames(Yk)[ty]]]= tidy(lm(scale(log2(Yk[,ty]+1))~mmp))
        #  glm_poisson(mmp,Yk[,ty])$be
     } 
    #------------------------------------------------------------------

       #per tissue joint hotspot betas
    hotspotBetasTis=rff[,hotspots$hotspot.markers]
    tissueHotspotBetas[[uc]]=hotspotBetasTis
    saveRDS(hotspotBetasTis, paste0(dout, 'tissue_joint_hotspotBetas.RDS'))

    #pool for permutations
    twgFDR= getFDRfx(rff, Yfe, Gf, cisMarkers, nperm=10)
    qtl.fdrs[[uc]]=twgFDR
    saveRDS(twgFDR, paste0(dout, 'FDRfx.RDS'))

    #per tissue cis
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(rff)),]
    cisBetasTis=rff[cbind(cMsubset$transcript, cMsubset$marker)]
    names(cisBetasTis)=cMsubset$transcript
    cisBetasTis=data.frame(transcript=names(cisBetasTis), sbeta=cisBetasTis,FDR=twgFDR[[4]](abs(cisBetasTis)), stringsAsFactors=F)
    cisBetasTis$FDR[is.na(cisBetasTis$FDR)]=1
    cisBetasTis$FDR[ cisBetasTis$FDR>1]=1
    tissueCisBetas[[uc]]=cisBetasTis
    saveRDS(cisBetasTis, paste0(dout, 'tissue_cisBetas.RDS'))

    
   
    # inputs 1) all peaks, 2) markers sorted into centimorgan bins, 3, fdr threshold
    tissue.hotspots = getHotspots(tcP, gmapd)
    saveRDS(tissue.hotspots,paste0(dout, 'within_tissue_hotspots.RDS'))
    withinTissueHotspots[[uc]]=tissue.hotspots
    
    if( length(tissue.hotspots[[1]])>0 ) {
        withinTissueHB=rff[,tissue.hotspots[[1]]]
        saveRDS(withinTissueHB,paste0(dout, 'within_tissue_hotspot_betas.RDS'))
        withinTissueHotspotBetas[[uc]]=withinTissueHB
    }
}




































#model total reads, batch, broad tissue classifications 
mm=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue) 
#+joint.covariates$fine_tissue)
BOLS=lm.fit(mm, log2(Yr+1))
Yresid=residuals(BOLS)
rownames(Yresid)=rownames(Yr)
colnames(Yresid)=colnames(Yr)
rm(BOLS)
#Yresid.var=colVars(Yresid,parallel=T)

Yre=standardise2(Yresid)
rownames(Yre)=rownames(Yresid)
colnames(Yre)=colnames(Yresid)

r=crossprod(Yre,G)/(nrow(Yre)-1)
#r2=crossprod(Yre, Gfit)/(nrow(Yre)-1)
saveRDS(r, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_std_beta.RDS')

LOD=-nrow(Yre)*log(1-r^2)/(2*log(10))
#saveRDS(r, '/data/single_cell_eQTL/elegans/results/joint_all_betas.RDS')
#saveRDS(LOD, '/data/single_cell_eQTL/elegans/results/joint_all_LODs.RDS')

#pool for permutations
wgFDR= getFDRfx(r, Yre, G, cisMarkers, nperm=10)
saveRDS(wgFDR, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_FDRfx.RDS')

cisBetasJoint=r[cbind(cisMarkers$transcript, cisMarkers$marker)]
names(cisBetasJoint)=cisMarkers$transcript

cisBetasJoint=data.frame(transcript=names(cisBetasJoint), 
                         sbeta=cisBetasJoint,
                         FDR=wgFDR[[4]](abs(cisBetasJoint)), stringsAsFactors=F)
cisBetasJoint$FDR[is.na(cisBetasJoint$FDR)]=1
cisBetasJoint$FDR[ cisBetasJoint$FDR>1]=1
saveRDS(cisBetasJoint, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_cisBetas.RDS')

cP=getGlobalPeaks(r, LOD,wgFDR,markerGR,transcript.data, chroms) 
saveRDS(cP, '/data/single_cell_eQTL/elegans/results/20191217/joint/jointPeaks.RDS') #XQTL_F4_2_jointPeaks_filt.RDS')

cP=readRDS('/data/single_cell_eQTL/elegans/results/20191217/joint/jointPeaks.RDS')


cP_nb=negbinRefit(cP, Yr,Gr, mm, fdr.thresh=.1)
saveRDS(cP_nb, '/data/single_cell_eQTL/elegans/results/20191217/joint/jointPeaks_withnegbin1.RDS')

rcp1=rbindlist(cP_nb,fill=T)
rcp1f=rcp1[!is.na(rcp1$negbin.p),]
rcp1f$FDR.negbin=p.adjust(rcp1f$negbin.p, method='fdr')

#cPf=rcp1[qvalue(rcp1$negbin.p)$qvalue<.01, ] #fthresh,]
#cPf=rcp1[(rcp1$negbin.p)<.001, ] #fthresh,]
#cPf=rcp1[p.adjust(rcp1$negbin.p, method='fdr')<.1, ] #fthresh,]

par(xaxs = "i", yaxs = "i")
#alpha=(-log10(FDR+1e-6)/6))
ggplot(cPf,aes(x=peakGpos,y=tGpos,alpha=-log10(negbin.p+1e-12)))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none') +
    scale_x_continuous(limits=c(0, max(gcoord.key)))+
    geom_hline(yintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=gcoord.key, color='lightblue', alpha=.9)+theme_classic()+
      theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

fthresh=.01
cPf=cP[cP$FDR<fthresh,]
ggplot(cPf,aes(x=peakGpos,y=tGpos,alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none') +
    scale_x_continuous(limits=c(0, max(gcoord.key)))+
    geom_hline(yintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=gcoord.key, color='lightblue', alpha=.9)+theme_classic()+
      theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())



# inputs 1) all peaks, 2) markers sorted into centimorgan bins, 3, fdr threshold
#hotspotsO = getHotspots(cP, gmapd, fdr.thresh=.05, pFDR=T)
hotspots = getHotspots(rcp1f, gmapd, fdr.thresh=.05, pFDR=F)

saveRDS(hotspots,'/data/single_cell_eQTL/elegans/results/20191217/joint/hotspots.RDS')
#hotspots=readRDS('/data/single_cell_eQTL/elegans/results/20191217/joint/hotspots.RDS')

hotspotBetasJoint=r[,hotspots$hotspot.markers]
saveRDS(hotspotBetasJoint, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_hotspotBetas.RDS')

######## Visualize --------------------------------------------------
fthresh=.01
cPf=cP[cP$FDR<fthresh,]
par(xaxs = "i", yaxs = "i")
ggplot(cPf,aes(x=peakGpos,y=tGpos,alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none') +
    scale_x_continuous(limits=c(0, max(gcoord.key)))+
    geom_hline(yintercept=gcoord.key, color='lightblue', alpha=.9)+
    geom_vline(xintercept=gcoord.key, color='lightblue', alpha=.9)+theme_classic()+
      theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
#---------------------------------------------------------------------

#[1] 133
#Error in x[, ii] : subscript out of bounds

#plot(markerGR$gcoord, af.tissue[[1]], xaxt='n', ylim=c(.2,.8), col=rainbow(20)[1], cex=.5, ylab='CB/N2+CB', xlab='')
#for(i in 2:length(af.tissue)){
#    points(markerGR$gcoord, af.tissue[[i]], xaxt='n', ylim=c(.2,.8), col=rainbow(20)[i], cex=.5)
#}
#legend('topright', legend=names(af.tissue),fill=rainbow(20)[1:20], cex=.75)



# Map within each tissue type --------------------------------------------------
qtl.fdrs=list()
tissueHotspotBetas=list()
tissueCisBetas=list()
tissueQTLPeaks=list()
tissueQTLPeaksNB=list()
withinTissueHotspots=list()
withinTissueHotspotBetas=list()

#qtl.maps=list()
#qtl.mlod=list()
#af.tissue=list()
#pc.lods=list()
#pc.betas=list()
#gets to be too much to retain these, should dump to files
#rho.matrices=list()
#residual.matrices=list()
fc=unique(as.character(joint.covariates$fine_tissue))
fc=fc[!is.na(fc)]
bc=unique(as.character(joint.covariates$broad_tissue))
bc=bc[!is.na(bc)]
#types=c(bc)
cty=c(rep('broad', length(bc)), rep('fine', length(fc)))
types=c(bc,fc)
#cty=c(rep('broad', length(bc)))

glms=list()
lms=list()
#error on 45
for(kk in 1:19) { 
    #20:86) { # 1:19 ) { #length(types)) {
    print(kk)
    uc=types[kk]
    print(uc)
    if(cty[kk]=='broad') {
        tk=which(joint.covariates$broad_tissue==uc)
        print('broad')
    }
    if(cty[kk]=='fine') {
        tk=which(joint.covariates$fine_tissue==uc & !is.na(joint.covariates$PC1))
        print('fine')
    }
    print(kk)
    if(length(tk)<6) { next; }
    covs=joint.covariates[tk,]
    Yk=Yr[tk,]
    #additional filter, expressed in at least 10 cells
    tcounts=Yk>0
    tcounts=colSums(tcounts)
    
    pick.filter=min(20, (length(tk)*.1))
    if(pick.filter<6) {pick.filter=6}
                    
    Yk=Yk[,tcounts>pick.filter]

    Gsub=Gr[tk,]
    #af.tissue[[uc]]=colSums(Gsub)/nrow(Gsub)
    
    Gf=standardise2(Gsub)
    colnames(Gf)=colnames(G)
    rownames(Gf)=rownames(Yk)

    #new poisson model ------------------------------------------
    mmp0=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch-1))
    mmp1=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch))

    #unstandardized genotypes 
    mmp=cbind(mmp0, Gsub[,hotspots$hotspot.markers])
       
    for(ty in 1:ncol(Yk)) {
       print(ty)
       glms[[uc]][[colnames(Yk)[ty]]]= tidy(glm(Yk[,ty]~mmp, family='poisson'))
       lms[[uc]][[colnames(Yk)[ty]]]= tidy(lm(scale(log2(Yk[,ty]+1))~mmp))
        #  glm_poisson(mmp,Yk[,ty])$be
     } 
    #------------------------------------------------------------------

    # OLS, replace with glm for next iteration
    # presumably run once only
    BOLS=lm(log2(Yk+1)~log2(covs$total)+covs$Batch)
    Yresid=residuals(BOLS)
    rownames(Yresid)=rownames(Yk)
    colnames(Yresid)=colnames(Yk)
    rm(BOLS)
    #residual.matrices[[uc]]=Yresid
    gc()
    #so slow ...
    Yfe=standardise2(Yresid)
    rownames(Yfe)=rownames(Yresid)
    colnames(Yfe)=colnames(Yresid)
       #transcripts that are singular after the regression (no variance = nothing to map on)
    shitT=unique(which(is.na(Yfe), arr.ind=T)[,2])
    if(length(shitT)>0) {    Yfe=Yfe[,-shitT] } 
  
    if(cty[kk]=='broad') {
        dout=paste0('/data/single_cell_eQTL/elegans/results/20191217/broad/', uc, '/')
    }
    if(cty[kk]=='fine') {
        dout=paste0('/data/single_cell_eQTL/elegans/results/20191217/fine/', uc, '/')
    }
    dir.create(dout)


    # then map within each classification
    rff=crossprod(Yfe,Gf)/(nrow(Yfe)-1)
    #saveRDS(rff, paste0(dout, 'tissue_std_beta.RDS'))

    #rho.matrices[[uc]]=rff
    LODr=-nrow(Yfe)*log(1-rff^2)/(2*log(10))
          
    #per tissue joint hotspot betas
    hotspotBetasTis=rff[,hotspots$hotspot.markers]
    tissueHotspotBetas[[uc]]=hotspotBetasTis
    saveRDS(hotspotBetasTis, paste0(dout, 'tissue_joint_hotspotBetas.RDS'))

    #pool for permutations
    twgFDR= getFDRfx(rff, Yfe, Gf, cisMarkers, nperm=10)
    qtl.fdrs[[uc]]=twgFDR
    saveRDS(twgFDR, paste0(dout, 'FDRfx.RDS'))

    #per tissue cis
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(rff)),]
    cisBetasTis=rff[cbind(cMsubset$transcript, cMsubset$marker)]
    names(cisBetasTis)=cMsubset$transcript
    cisBetasTis=data.frame(transcript=names(cisBetasTis), sbeta=cisBetasTis,FDR=twgFDR[[4]](abs(cisBetasTis)), stringsAsFactors=F)
    cisBetasTis$FDR[is.na(cisBetasTis$FDR)]=1
    cisBetasTis$FDR[ cisBetasTis$FDR>1]=1
    tissueCisBetas[[uc]]=cisBetasTis
    saveRDS(cisBetasTis, paste0(dout, 'tissue_cisBetas.RDS'))

    tcP=getGlobalPeaks(rff, LODr,twgFDR, markerGR,transcript.data, chroms) 
    tissueQTLPeaks[[uc]]=tcP
    saveRDS(tcP, paste0(dout,'Peaks.RDS')) #XQTL_F4_2_jointPeaks_filt.RDS')

    tcP_nb=negbinRefit(tcP, Yk,Gsub, mmp1, fdr.thresh=.1)

    trcp1=rbindlist(tcP_nb,fill=T)
    trcp1f=trcp1[!is.na(trcp1$negbin.p),]
    trcp1f$FDR.negbin=p.adjust(trcp1f$negbin.p, method='fdr')

   
    # inputs 1) all peaks, 2) markers sorted into centimorgan bins, 3, fdr threshold
    tissue.hotspots = getHotspots(tcP, gmapd)
    saveRDS(tissue.hotspots,paste0(dout, 'within_tissue_hotspots.RDS'))
    withinTissueHotspots[[uc]]=tissue.hotspots
    
    if( length(tissue.hotspots[[1]])>0 ) {
        withinTissueHB=rff[,tissue.hotspots[[1]]]
        saveRDS(withinTissueHB,paste0(dout, 'within_tissue_hotspot_betas.RDS'))
        withinTissueHotspotBetas[[uc]]=withinTissueHB
    }
}
perTissueResults=list(
qtl.fdrs                =qtl.fdrs,
tissueHotspotBetas       =tissueHotspotBetas  ,     
tissueCisBetas           = tissueCisBetas      ,    
tissueQTLPeaks          =  tissueQTLPeaks       ,   
withinTissueHotspots    =  withinTissueHotspots  ,  
withinTissueHotspotBetas =  withinTissueHotspotBetas)
saveRDS(perTissueResults, file='/data/single_cell_eQTL/elegans/results/20191217/perTissueLists.RDS')

perTissueResults=readRDS('/data/single_cell_eQTL/elegans/results/20191217/perTissueLists.RDS')

ptrq=perTissueResults$tissueQTLPeak
ptrq=ptrq$"Body Wall Muscle"
ptrq=ptrq[ptrq$FDR<.1,]
ptrq=split(ptrq,ptrq$transcript)
ptrq=ptrq[sapply(ptrq,nrow)>0]


gxy=lapply(glms, function(x) sapply(x, function(y) y$estimate))
lxy=lapply(lms, function(x) sapply(x, function(y) y$estimate))
saveRDS(gxy, file = '/data/single_cell_eQTL/elegans/results/20191217/hotspot_betas_glm_pois.RDS')
saveRDS(lxy, file = '/data/single_cell_eQTL/elegans/results/20191217/hotspot_betas_lm.RDS')



lbeta=lxy[['Body Wall Muscle']][,7:17]
gbeta=gxy[['Body Wall Muscle']][,7:17]

lbeta=lxy[['Germline']][,7:17]
gbeta=gxy[['Germline']][,7:17]







qtl.fdrs=list()
tissueHotspotBetas                =list()
tissueCisBetas                    =list()
tissueQTLPeaks                    =list()
withinTissueHotspots              =list()
withinTissueHotspotBetas          =list()



qm2=rbindlist(tissueQTLPeaks,idcol='tissue')

#qm$peakGpos=markerGR$gcoord[match(qm$peak.marker, colnames(G))]
#qm$tGpos=transcript.data$gcoord[match(qm$transcript, transcript.data$wormbase_gene)]
#qm$peakLGpos=markerGR$gcoord[match(qm$CI.l, colnames(G))]
#qm$peakRGpos=markerGR$gcoord[match(qm$CI.r, colnames(G))]

#qm$markerpos=match(qm$peak.marker, colnames(G))
#qm$tpos=match(qm$transcript, transcript.data$wormbase_gene)
tkk=names(which(sapply(tissueQTLPeaks, function(x) sum(x$FDR<.05))>50))
qm2=qm2[qm2$tissue %in% tkk,]
qm2=qm2[qm2$FDR<.01,]
ggplot(qm2,aes(x=peakGpos,y=tGpos, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    facet_wrap(~tissue)+geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

    #}, error=function(e) {NULL;} )

 ##map on pcs
    #if(cty[kk]=='broad' & uc!='Neuron') {
    #    pc.norm1=scale(residuals(lm(covs$PC1~log2(covs$total)+covs$Batch)))
    #    pc.norm2=scale(residuals(lm(covs$PC2~log2(covs$total)+covs$Batch)))
    #    Ypc=cbind(pc.norm1, pc.norm2)
    #    colnames(Ypc)=c('PC1', 'PC2')

    #     rff=crossprod(Ypc,Gf)/(nrow(Ypc)-1)
    #     #rho.matrices[[uc]]=rff
    #     LODr=-nrow(Ypc)*log(1-rff^2)/(2*log(10))
    #     pc.lods[[uc]]=LODr
    #     pc.betas[[uc]]=rff
    #}

    #if(cty[kk]=='fine' ) {
    #    pc.norm1=scale(residuals(lm(covs$PC1~log2(covs$total)+covs$Batch)))
    #    #pc.norm2=scale(residuals(lm(covs$PC2~log2(covs$total)+covs$Batch)))
    #    Ypc=cbind(pc.norm1)#, pc.norm2)
    #    colnames(Ypc)=c('pseudotime') #PC1', 'PC2')
    #     rff=crossprod(Ypc,Gf)/(nrow(Ypc)-1)
    #     #rho.matrices[[uc]]=rff
    #     LODr=-nrow(Ypc)*log(1-rff^2)/(2*log(10))
    #     pc.lods[[uc]]=LODr
    #     pc.betas[[uc]]=rff
    #}
 #  #multivariate scan for hotspots 
 #   mLOD=list()
 #   for(cc in chroms) {
 #       print(cc)
 #       ng=transcript.data$wormbase_gene[transcript.data$chromosome_name!=cc]
 #       YreS=Yfe[,(colnames(Yfe) %in% ng)]
 #       moi=which(as.character(seqnames(markerGR))==cc)
 #       #test=hd.eigen(t(Yre), center=F, scale=F,vectors=T)
 #       pc.to.retain=40
 #       yreduced=hd.eigen(t(YreS), center=F, scale=F, k=pc.to.retain, vectors=T)
 #       testF=mvn.scanone(Gf[,moi], yreduced$vectors)
 #       testN=determinant(crossprod(yreduced$vectors), logarithm=T)$modulus
 #       mLODv=(nrow(yreduced$vectors)/2)*(testN-testF)/(2*log(10))
 #       mLOD[[cc]]=mLODv
 #   }
 #   mLOD=do.call('c', mLOD)    
 #   qtl.mlod[[uc]]=mLOD

 #   plot(markerGR$gcoord, mLOD, main=uc)
 #   abline(v=gcoord.key, lty=2, color='lightblue')





# cis model 
#mm=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue) #+joint.covariates$fine_tissue)
#BOLS=lm.fit(mm, log2(Yr+1), intercept=F)
library(broom)

# for each transcript, count up the number of cells per tissue the transcript is expressed in with at least one read
pte=apply(Yr, 2, function(x) sapply(split(x, joint.covariates$broad_tissue), function(y) sum(y>0)))
# hard cutoff for 'is expressed in' each broad tissue
t.cutoff=round(table(joint.covariates$broad_tissue)*.1)
t.cutoff20=20
is.expressed10=apply(pte,2,function(x) x>t.cutoff)
is.expressed20=apply(pte,2,function(x) x>t.cutoff20)

saveRDS(is.expressed10, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_10_percent_of_broad_tissue.RDS')
saveRDS(is.expressed20, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_20_cells_broad_tissue.RDS')

# build model matrix for broad tissue assignment
btm=model.matrix(log2(Yr[,1]+1)~-1+broad_tissue, data=joint.covariates)
colnames(btm)=gsub('broad_tissue', '', colnames(btm))

#cResid=list()
#c.model=list()
# precompute some permutation vectors -----------------
set.seed(30)
permlist=list()
for(perm in 1:10){
   print(perm)
   pchar=as.character(perm)
   pvec=sample(1:nrow(Yr))
   permlist[[pchar]]=pvec
}
permlist=do.call('cbind', permlist)
#------------------------------------------------------
#library(DescTools)


# extract matrix of hotspot markers 
hm=G[,hotspots$hotspot.markers]
# build matrix of hotspot-tissue interaction covariates 
hml=list()
for(i in 1:ncol(hm)){
    hml[[colnames(hm)[i]]]=hm[,i]*btm
    colnames(hml[[colnames(hm)[i]]])=paste0(colnames(hm)[i], ':', colnames(btm))
}
hml=scale(do.call('cbind', hml))

# build design matrix for counts, batch, and broad_tissue 
m0=model.matrix(lm(log2(Yr[,1]+1)~log2(total)+Batch+broad_tissue,data=joint.covariates))
# n - 1 for calculating pearson R
nf=nrow(hml)-1

# flag all transcripts classified as expressed in at least one tissue 
#expressedTranscripts=names(which(colSums(is.expressed10)>0))
expressedTranscripts=names(which(colSums(is.expressed20)>0))

# filter cisMarkers based on transcript expressed in at least one tissue 
cisMarkersFiltered=cisMarkers[cisMarkers$transcript %in% expressedTranscripts,]

#is.expressedFiltered=is.expressed10[,colnames(is.expressed10) %in% expressedTranscripts]
is.expressedFiltered=is.expressed20[,colnames(is.expressed20) %in% expressedTranscripts]

# lists to hold results of interactions tests 
c.int.r=list()
h.int.r=list()

# just test for interaction effects
for(gi in 1:nrow(cisMarkersFiltered)){
        print(gi)

        # get transcript and cis marker
        transcriptY=cisMarkersFiltered[gi,1]
        cMarker=cisMarkersFiltered[gi,2]
        cMarkerG=G[,cMarker]
        # append marginal cis effet to existing design matrix
        XX=cbind(m0,cMarkerG)
       
        m1=.lm.fit(XX,log2(Yr[,transcriptY]+1)) 
        r1=scale(residuals(m1))
        rp1=apply(permlist, 2, function(x) r1[x])  
        
        # construct covariate for cis-tissue interaction covariates
        bm=scale(btm*cMarkerG)
        # test
        c1=crossprod(cbind(r1,rp1),bm)/nf
        c.int.r[[transcriptY]]=c1 
  
        # add cis-interaction terms to model
        XX2=cbind(XX,bm)
        m2=.lm.fit(cbind(XX2,hm),log2(Yr[,transcriptY]+1)) 
        r2=scale(residuals(m2))
        rp2=apply(permlist, 2, function(x) r2[x])  
        
        # test for hotspot-tissue interactions 
        c4=crossprod(cbind(r2,rp2),hml)/nf
        h.int.r[[transcriptY]]=c4        
}
c.int.obs=sapply(c.int.r, function(x) x[1,])
c.int.exp=abind(lapply(c.int.r, function(x) x[2:11,]), along=3)

c.int.obs.stat=abs(c.int.obs)
c.int.exp.stat=abs(c.int.exp)

vint=seq(0.001,max(c.int.obs.stat)+.001, .001)
ctisFDR=c.int.obs
ctisFDR[is.numeric(ctisFDR)]=NA

for(tis in rownames(c.int.obs) ){
    print(tis)
    tstat=c.int.obs.stat[tis,]
    #names(is.expressedFiltered[tis,])
    ie=(names(tstat)%in%names(which(is.expressedFiltered[tis,])))
    tisexp=tstat[ie]

    obsPcnt=sapply(vint, function(thresh) { sum(tisexp>thresh) } ) #apply(tisexp,1, function(x) sum(x>thresh) ) } )
    names(obsPcnt)=vint

    expPcnt=colMeans(sapply(vint, function(thresh) { 
               apply(c.int.exp.stat[,tis,ie],1, function(x) sum(x>thresh))
           }))
    names(expPcnt)=vint
    
   pFDR = expPcnt/obsPcnt
   pFDR[is.na(pFDR)]=0
   pFDR[!is.finite(pFDR)]=0
   #to make sure this is monotonic
   pFDR=rev(cummax(rev(pFDR)))

   pFDR[1]=pFDR[1]-(1e-3)
   fdrFX=approxfun(pFDR, vint, ties=mean) 
   rtoFDRfx=approxfun(vint,pFDR, ties=mean) 
   tfdr=rtoFDRfx(tisexp)
   tfdr[is.na(tfdr)]=1
   ctisFDR[tis,ie]=tfdr
}

cint=list(cisTissueBetas=c.int.obs, cisTissueFDR=ctisFDR)
#saveRDS(cint, '/data/single_cell_eQTL/elegans/results/20191217/tissue_interactions/cis_tissue_interaction_betas_10.RDS')
saveRDS(cint, '/data/single_cell_eQTL/elegans/results/20191217/tissue_interactions/cis_tissue_interaction_betas_20.RDS')


h.int.obs=sapply(h.int.r, function(x) x[1,])
h.int.exp=abind(lapply(h.int.r, function(x) x[2:11,]), along=3)

h.int.obs.stat=abs(h.int.obs)
h.int.exp.stat=abs(h.int.exp)

htisFDR=h.int.obs
htisFDR[is.numeric(htisFDR)]=NA

for(tis in rownames(h.int.obs) ){
    print(tis)
    tstat=h.int.obs.stat[tis,]
    toi=strsplit(tis, ':')[[1]][2]
    #names(is.expressedFiltered[tis,])
    ie=(names(tstat)%in%names(which(is.expressedFiltered[toi,])))
    tisexp=tstat[ie]

    obsPcnt=sapply(vint, function(thresh) { sum(tisexp>thresh) } ) #apply(tisexp,1, function(x) sum(x>thresh) ) } )
    names(obsPcnt)=vint

    expPcnt=colMeans(sapply(vint, function(thresh) { 
               apply(h.int.exp.stat[,tis,ie],1, function(x) sum(x>thresh))
           }))
    names(expPcnt)=vint
    
   pFDR = expPcnt/obsPcnt
   pFDR[is.na(pFDR)]=0
   pFDR[!is.finite(pFDR)]=0
   #to make sure this is monotonic
   pFDR=rev(cummax(rev(pFDR)))

   pFDR[1]=pFDR[1]-(1e-3)
   fdrFX=approxfun(pFDR, vint, ties=mean) 
   rtoFDRfx=approxfun(vint,pFDR, ties=mean) 
   tfdr=rtoFDRfx(tisexp)
   tfdr[is.na(tfdr)]=1
   htisFDR[tis,ie]=tfdr
}
hint=list(hotTissueBetas=h.int.obs, hotTissueFDR=htisFDR)
#saveRDS(hint, '/data/single_cell_eQTL/elegans/results/20191217/tissue_interactions/hotspot_tissue_interaction_betas_10.RDS')
saveRDS(hint, '/data/single_cell_eQTL/elegans/results/20191217/tissue_interactions/hotspot_tissue_interaction_betas_20.RDS')


hsig=apply(htisFDR,1,function(x) sum(x<.05,na.rm=T))


dhsig= data.frame(tstrsplit(names(hsig),':'), cnt=as.vector(hsig))
names(dhsig)[c(1,2)]=c('hotspot', 'tissue')

ggplot(dhsig, aes(x=hotspot, y=cnt))+geom_bar(stat='identity')+facet_wrap(~tissue)

#c.int.obs[tis,is.na(rtoFDRfx[[tis]](abs(c.int.obs[tis,])))]




























t.int.Z=sapply(c.int.model, function(x) x[,'stat'])
is.expressed2=is.expressed
is.expressed2=is.expressed2[,colnames(t.int.Z)] #t

t.int.Z[!is.expressed2]=0
t.int.P=sapply(c.int.model, function(x) x[,'pvalue'])
t.int.P[!is.expressed2]=NA

t.int.Z=t.int.Z

c.model.perm=list()
c.int.model.perm=list()


mLODj=list()
#joint scan for hotspots
for(cc in chroms[-1]) {
    print(cc)
    cPsigm=cPsig[cPsig$tchr!=cc,]

    #ng=transcript.data[transcript.data$chromosome_name !='I',]
    YreS=Yre[,colnames(Yre) %in% unique(cPsigm$transcript)] #ng$wormbase_gene]

    moi=grep(paste0('^',cc,'_'), colnames(G))
    pc.to.retain=40
    yreduced=hd.eigen(t(YreS), center=F, scale=F, k=pc.to.retain, vectors=T)
    
    testF=mvn.scanone(G[,moi], yreduced$vectors)
    testN=determinant(crossprod(yreduced$vectors), logarithm=T)$modulus
    mLODv=(nrow(yreduced$vectors)/2)*(testN-testF)/(2*log(10))
    mLODj[[cc]]=mLODv

}



#ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/100919/elegans_joint_CI.png')




#cPt=do.call('rbind',cPeaksT)

#plot(match(cPt$peak.marker[cPt$FDR<.9], colnames(G)),
#     match(cPt$transcript[cPt$FDR<.9], transcript.data$wormbase_gene), 
#      xlab='marker index', ylab='transcript index', main='joint analysis FDR < 1%')




#ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/100919/elegans_tissue_split.png')


#+scale_color_brewer(palette='Set1')


ggplot(qm,aes(x=peakGpos,y=tGpos,color=tissue, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x = peakLGpos, y = tGpos, xend = peakRGpos, yend = tGpos)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/100919/elegans_tissue_rainbow.png')













































#new genotypes 
r2=crossprod(Yre,Gr1b)/(nrow(Yre)-1)


cis1=-nrow(Yre)*log(1-r1[cis.index]^2)/(2*log(10))
cis2=-nrow(Yre)*log(1-r2[cis.index]^2)/(2*log(10))



cisP_G=rep(NA, nrow(cisMarkers))
cisP_G2=rep(NA, nrow(cisMarkers))
for(rr in 1:nrow(cisMarkers)) {
    print(rr)
    lmm=lm(Yre[,cisMarkers$transcript[rr]]~G2[crit,cisMarkers$marker[rr]])
    plot(G2[crit,cisMarkers$marker[rr]], Yre[,cisMarkers$transcript[rr]],col="#00000022") 
    abline(lm(Yre[,cisMarkers$transcript[rr]]~G2[crit,cisMarkers$marker[rr]]), col='red')
    plot(G2[crit,cisMarkers$marker[rr]], Y[crit,cisMarkers$transcript[rr]],col="#00000022") 

    cisMarker=G2[crit,cisMarkers$marker[rr]]
    cisTranscriptL=log2(Y[crit,cisMarkers$transcript[rr]]+1)
    cisTranscriptR=(Y[crit,cisMarkers$transcript[rr]])

    cisTissues=joint.covariates$broad_tissue
    lmm=lm(cisTranscriptL~log2(joint.covariates$total)+joint.covariates$Batch+cisMarker*cisTissues)

    glmm=glm.nb(cisTranscriptR~offset(log(joint.covariates$total))+joint.covariates$Batch+cisMarker*cisTissues)
    
    
    G2[crit,cisMarkers$marker[rr]])


    cisP_G2[rr]=(drop1(lmm, test='Chisq'))[[5]][2]
    lmm=lm(Yre[,cisMarkers$transcript[rr]]~rank(G[crit,cisMarkers$marker[rr]]) )
    cisP_G[rr]=(drop1(lmm, test='Chisq'))[[5]][2]
}












# why is this so fucking slow 
Gr1=standardise2(colRanks(G[crit,],parallel=T))
colnames(Gr1)=colnames(G)
rownames(Gr1)=rownames(Yre)

r=crossprod(Yre,Gr1)/(nrow(Yre)-1)
LOD1=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL1=which(LOD1>8, arr.ind=T)
rm(r)
rm(LOD1)
plot(sigL1[,2], sigL1[,1], pch=21,main='ranks r/qtl genotypes')
#------------------

#-------------------------------------------------
Gr2=standardise2(colRanks(G2[crit,],parallel=T))
colnames(Gr2)=colnames(G)
rownames(Gr2)=rownames(Yre)

r=crossprod(Yre,Gr2)/(nrow(Yre)-1)
LOD2=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL2=which(LOD2>8, arr.ind=T)
plot(sigL2[,2], sigL2[,1], pch=21,main='ranks custom genotypes')

rm(r)
rm(LOD2)
#-------------------------------------------------


#-------------------------------------------------
Gr2b=standardise2(G2[crit,])
colnames(Gr2b)=colnames(G)
rownames(Gr2b)=rownames(Yre)
r=crossprod(Yre,Gr2b)/(nrow(Yre)-1)
LOD3=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL3=which(LOD3>8, arr.ind=T)
plot(sigL3[,2], sigL3[,1], pch=21,main='custom genotypes')

rm(r)
rm(LOD3)
#---------------------------------------------------

#------------------------------------------
Gr1b=standardise2(G[crit,])
colnames(Gr1b)=colnames(G)
rownames(Gr1b)=rownames(Yre)
r=crossprod(Yre,Gr1b)/(nrow(Yre)-1)
LOD4=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL4=which(LOD4>8, arr.ind=T)
rm(r)
rm(LOD4)
x11()
plot(sigL4[,2], sigL4[,1], pch=21,main='r/qtl genotypes')


par(mfrow=c(2,2))
plot(sigL1[,2], sigL1[,1], pch=21, cex=.1, col="#00000022", main='ranks r/qtl genotypes',xlab='marker index', ylab='transcript index')
plot(sigL4[,2], sigL4[,1], pch=21,, cex=.1, col="#00000022", main='r/qtl genotypes',xlab='marker index', ylab='transcript index')
plot(sigL2[,2], sigL2[,1], pch=21,, cex=.1, col="#00000022", main='ranks custom genotypes',xlab='marker index', ylab='transcript index')
plot(sigL3[,2], sigL3[,1], pch=21,, cex=.1, col="#00000022", main='custom genotypes',xlab='marker index', ylab='transcript index')



G2sum=colSums(G2)
G2rs=rowSums(G2)



tt=eMats/colSums(eMats)
test=eLikeHMM(tt,ttMat)


#CB=N2
#n=5
#k=1
#nmk=n-k
#bco=choose(n,k)
#ee=.0001

#bco*((1-ee)^k)*(ee^(nmk))
#bco*(.5)^n
#bco*((1-ee)^nmk)*ee^k

#ref/het/alt
#pRR  = dbinom(n-k,n,ee) *      (1-r)/2
#pHet = choose(n,k)*(1/(2^n))*   r
#pAA  = dbinom(k,n,ee) *        (1-r)/2


nind=nrow(cross2$cross_info)
cc=cut(1:nind, 100)

f='1'
gps=readRDS(paste0('/data/single_cell_eQTL/elegans/XQTL_F4_2/hmm_v2/', f))

#501 is weird

for(h in 600:691){
    i=rownames(gps[[1]])[h]# 691
    coii='II'
    coi=paste0('^', coii, '_') #'^II_'

    #plot(gps[[2]][i,1,])
    N_bs=N2.counts[grep(coi, rownames(N2.counts)),i]
    C_bs=CB.counts[grep(coi, rownames(CB.counts)),i]

   # N_bs0=N2.counts0[grep(coi, rownames(N2.counts0)),i]
   # C_bs0=CB.counts0[grep(coi, rownames(CB.counts0)),i]

   # N_bsR=N2.countsR[grep(coi, rownames(N2.countsR)),i]
   # C_bsR=CB.countsR[grep(coi, rownames(CB.countsR)),i]

    cGcall=cross2$geno[[coii]][i,]

    par(mfrow=c(2,1))
    plot(.5*gps[[coii]][i,2,]+gps[[coii]][i,3,],main='p(C)', col='blue', type='l')
    #points(gps[[coii]][i,2,],main='p(NC)', col='grey')
    #points(gps[[coii]][i,3,],main='p(CC)', col='blue')
    #test=unlist(cGcall, function(x) switch(x, 'blue', 'grey', 'orange', 'purple', 'green'))
    cme=which(N_bs>0)
    if(length(cme)==0) {next;}

    plot(cme,N_bs[cme], type='h', ylim=c(-10,10), main='br, strict',col=cGcall[names(cme)]) #sapply(cGcall, function(x) switch(x, 'blue', 'grey', 'orange', 'purple', 'green')))
    points(cme, -C_bs[cme], type='h', col=cGcall[names(cme)]) #, col=cGcall)
    abline(h=0)

   # plot(-N_bs0, type='h', ylim=c(-10,10), main='br')
   # points(C_bs0, type='h') #, col=cGcall)

    #plot(-N_bsR, type='h', ylim=c(-10,10), main='raw')
    #points(C_bsR, type='h') #, col=cGcall)

    print(h)
    readline()
    }
#black is CC
#red is NC
#green is NN
#blue is NN or NC
#cyan is NC or CC





hets=which(cGcall==2)
points(hets, rep(0,length(hets)), col='red', pch=17, cex=2)
notNN=which(cGcall==5)
points(notNN, rep(0,length(notNN)), col='cyan',pch=17,cex=2)
notCC=which(cGcall==4)
points(notCC, rep(0,length(notCC)), col='blue', pch=17, cex=2)
isCC=which(cGcall==1)
points(isCC, rep(0,length(isCC)), col='black', pch=17, cex=2)
isNN=which(cGcall==3)
points(isNN, rep(0,length(isNN)), col='green',pch=17, cex=2)












#seq_along(covariates$broad_tissue), paste(covariates$Batch covariates$broad_tissue)
set.seed(30)

df1=data.frame(tis=joint.covariates$broad_tissue, batch=joint.covariates$Batch, id=seq_along(joint.covariates$broad_tissue))
df1s=split(df1, paste(df1$tis,df1$batch ))
df1s=lapply(df1s, function(x){ x$pid=sample(x$id); return(x)})
pvec=do.call('c', lapply(df1s, function(x) x$pid))

for(perm in 1:5) {
    #perm=1
    m0=lm(log2(Yr[,1]+1)~log2(total)+Batch+broad_tissue,data=joint.covariates)
    X0=model.matrix(m0)
    #bm=bm[pvec,]
   
    #pvec=as.vector(unlist(sapply(split(seq_along(covariates$broad_tissue), covariates$broad_tissue),sample)))
    for(gi in 1:nrow(cisMarkers)){
        print(gi)
        transcriptY=cisMarkers[gi,1]
        cMarker=cisMarkers[gi,2]
        cMarkerG=G[,cMarker]
        XX=cbind(X0,cMarkerG)
        XX=XX[pvec,]
        #XX=XX[pvec,]
        XX[,'log2(total)']=X0[,'log2(total)']
        bm=btm*cMarkerG

        yp=log2(Yr[,transcriptY]+1)

        m1=.lm.fit(XX,yp) #log2(Yr[,transcriptY]+1)) 
        m3=lm(yp~XX) #log2(total)+Batch+broad_tissue+cMarkerG,data=covariates)

        r1=residuals(m1)
        t2=univglms(r1, bm[pvec,])


        m0=lm(log2(Yr[,transcriptY]+1)~log2(total)+Batch+broad_tissue+cisMarker,data=covp) #ariates)
        #r0=residuals(m0)
        #cResid[[transcriptY]]=r0
        c.model.perm[[as.character(perm)]][[transcriptY]]=tidy(m0)
        # this is slow
        #cResid[,transcriptY]=r0 
        #esiduals(m0)
        m1=update(m0,~.+broad_tissue:cisMarker, data=covariates)
        c.int.model.perm[[as.character(perm)]][[transcriptY]]=tidy(m1)
    }
}








    #m2=lm(r0~broad_tissue:cisMarker,data=covariates)

    #pvec=as.vector(unlist(sapply(split(seq_along(covariates$broad_tissue), covariates$broad_tissue),sample)))
    #covp=covariates[pvec,]
    #mp=lm(log2(Yr[pvec,transcriptY]+1)~log2(total)+Batch+broad_tissue:cisMarker,data=covp)


    #m00=glm.nb(Yr[,transcriptY]~log2(total)+Batch+broad_tissue+cisMarker,data=covariates)
    #m11=update(m00,~.+broad_tissue:cisMarker, data=covariates) 
    #glm.nb(Yr[,transcriptY]~log2(total)+Batch+broad_tissue+cisMarker,data=covariates)



summary(glm.nb(Yr[,transcriptY]~offset(log(total))+Batch+broad_tissue*cisMarker, data=covariates))
summary(lm(log2(Yr[,transcriptY]+1)~log2(total)+Batch+broad_tissue*cisMarker,data=covariates)) #+broad_tissue*hmarkers, data=covariates))



summary(lm(log2(Yr[,transcriptY]+1)~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue*cMarkerG))
summary(lm(log2(Yr[,transcriptY]+1)~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue*cMarkerG))

summary(glm.nb(Yr[,transcriptY]~log(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue*cMarkerG))
summary(glm.nb(Yr[,transcriptY]~offset(log(total))+Batch+broad_tissue*cisMarker+broad_tissue*I_184284, data=test))



























    fdrfx.fc=getFDRfx(rff,Yfe,Gf,nperm=5) #,vint=seq(0.001,.45,.001))
    
    
    qtl.fdrs[[uc]]=fdrfx.fc
   
    cPeaksT=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(rff))
        rS=rff[,moi]
        #mstat=apply(rS,1,max)
        mstat.pos=apply(abs(rS),1,which.max)
        lookup=cbind(1:length(mstat.pos), mstat.pos)
        mstat=rS[lookup]
        LODstat=LODr[,moi][lookup]
        CIs=matrix(NA,length(mstat),2)
        cmoi=colnames(LODr[,moi])
        #this could use a speedup
        for(peak in 1:length(mstat)){
            LODrv=LODr[peak,moi]
            CIs[peak,]=cmoi[range(which(LODrv>LODstat[peak]-1.5))]
        }   
        cPeaksT[[cc]]=data.frame( transcript=rownames(rS),
                             peak.marker=colnames(rS)[mstat.pos],
                             CI.l=CIs[,1],
                             CI.r=CIs[,2],
                             LOD=LODstat, 
                             sbeta=mstat,
                             FDR=fdrfx.fc[[2]](abs(mstat)) )
        cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
        cPeaksT[[cc]]$FDR[(cPeaksT[[cc]]$FDR)>1]=1
    }
 
   qtl.maps[[uc]]=rbindlist(cPeaksT, idcol='chrom')

   cPt=qtl.maps[[uc]]
   if(sum(cPt$FDR<.2)>0) {
       plot(match(cPt$peak.marker[cPt$FDR<.2], colnames(G)),
         match(cPt$transcript[cPt$FDR<.2], transcript.data$wormbase_gene), 
          xlab='marker index', ylab='transcript index', main=paste(uc, 'joint analysis FDR < 20%'))
   }
   print(sum(qtl.maps[[uc]]$FDR<.2, na.rm=T))
}
#---------------------------------------------------------------------------------------------------

#saveRDS(qtl.maps ,  file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePeaks_filt.RDS'               )
#saveRDS(pc.lods  ,  file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePCs_PseudoTime_filt.RDS'      )
#saveRDS(qtl.mlod ,  file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissueMV_Peaks_filt.RDS'            )
#saveRDS(pc.betas ,  file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePCs_PseudoTime_Betas_filt.RDS')

qtl.maps=readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePeaks_filt.RDS'               )
pc.lods =readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePCs_PseudoTime_filt.RDS'      )
#qtl.mlod=readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissueMV_Peaks_filt.RDS'            )
pc.betas=readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_perTissuePCs_PseudoTime_Betas_filt.RDS')


plot(log2(colSums(Ns)), log2(colsSums(Cs)), xlab='log2(informative N2 sites)

plot(log2(rowSums(Ns>0)), ylim=c(-5,5), type='h')
points(-log2(rowSums(Cs>0)), type='h')
abline(v=sapply(split(1:nrow(Ns),gmapd[gmapd$chrom=='X',]$cm5bin), max),col='red', lwd=3)

split(rowSums(Ns>0), gmapd[gmapd$chrom=='X',]$cm5bin)

