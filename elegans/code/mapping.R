install.packages('pracma')
install.packages('Rfast2')
install.packages('fastglm')
install.packages('ggsci')
install.packages('caret')
#BiocManager::install("qvalue")
#install.packages('float')

library(data.table)
library(ggsci)
library(ggplot2)
library(pracma)
library(Rfast2)
library(Rfast)
library(fastglm)
library(MASS)
library(caret)
#library(float)
source('/data/single_cell_eQTL/elegans/code/utilities.R')
source('/data/single_cell_eQTL/elegans/code/mapping_fxs.R')

# cells with NAs for broad_tissue classification are useless here
#crit=!is.na(fine.classification$broad) #& !is.na(fine.classification$fine) 
# & log2(fine.classification$total)>8.5 & log2(fine.classification$total)<15.5
crit=!is.na(joint.covariates$broad_tissue)
joint.covariates=joint.covariates[crit,]
saveRDS(joint.covariates, '/data/single_cell_eQTL/elegans/results/20191217/covariates.RDS')

#filter count matrix ------------------------------------------------------
Y=t(as.matrix(counts))
Yr=Y[crit,]

#calc Variance for each transcript, if zero, boot it
Yr.var=colVars(Yr, parallel=T)
Yr=Yr[,Yr.var>0]
#additional filter, expressed in at least 21 cells
tcounts=Yr>0
tcounts=colSums(tcounts)
Yr=Yr[,tcounts>20]
saveRDS(Yr, '/data/single_cell_eQTL/elegans/results/20191217/transcript_matrix.RDS')
# -----------------------------------------------------------------------------------------

# filter genotypes ------------------------------------------------------------------------
#geno=do.call('cbind', geno)
Gr=geno[crit,]
saveRDS(Gr, '/data/single_cell_eQTL/elegans/results/20191217/geno_matrix.RDS')
G=standardise2(Gr)

grs=split(data.frame(Gr), joint.covariates$broad_tissue)
gvar=lapply(grs, function(x) colVars(data.matrix(x), parallel=T))
saveRDS(gvar, '/data/single_cell_eQTL/elegans/results/20191217/broad/MarkerVariance.RDS')
#-----------------------------------------------------------------------------------------

#correlate genotypes across strains 
#A=tcrossprod(G)
#A=A/(nrow(G)-1)
#rr=which(A>.99, arr.ind=T)
#rr=rr[rr[,1]!=rr[,2],]

# bin the genome into 5cm bins -------------------------------------------------------------------------
gmapd=binMarkersToGeneticMap(Gr,gmap, bin.size=10)
#saveRDS(gmapd, file = '/data/single_cell_eQTL/elegans/results/20191217/markers_to_5cm_bins.RDS')
#-------------------------------------------------------------------------------------------------------

# identify closest marker to transcript for cis only test (to do) ----------------------------------
g2m=tstrsplit(colnames(G),'_')
# would functionalize this, but using markerGR later
markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*",name=colnames(G)))
markerGR$gcoord=gcoord.key[as.character(seqnames(markerGR))]+start(markerGR)

#tag=transcript.annotations[transcript.annotations$type=='gene',]
#tag=tag[match(colnames(Yr), tag$gene_id),]
#cisMarkers=data.frame(transcript=colnames(Yr), marker=colnames(G)[nearest(tag,markerGR)],stringsAsFactors=F)
#cisMarkers=na.omit(cisMarkers)
#saveRDS(cisMarkers, '/data/single_cell_eQTL/elegans/results/20191217/closestMarker.RDS')
#-------------------------------------------------------------------------------------------------

#per tissue betas  (no QTL) -----------------------------------------------------------------------
#m0=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch)
#r0=(residuals(lm.fit(m0,log2(Yr+1))))
#m1=model.matrix(Yr[,1]~joint.covariates$broad_tissue-1)
#ctm=lm.fit(m1, r0)
#ctm=coef(ctm)
#rownames(ctm)=gsub('joint.covariates\\$broad_tissue', '', rownames(ctm))
#saveRDS(ctm, file = '/data/single_cell_eQTL/elegans/results/20191217/transcript_tissue_betas.RDS')
#---------------------------------------------------------------------------------------------------


#if interrupted before line 82----------------------------------------------------------------
joint.covariates=readRDS( '/data/single_cell_eQTL/elegans/results/20191217/covariates.RDS')
Yr=readRDS('/data/single_cell_eQTL/elegans/results/20191217/transcript_matrix.RDS')
Gr=readRDS('/data/single_cell_eQTL/elegans/results/20191217/geno_matrix.RDS')
gmapd=readRDS('/data/single_cell_eQTL/elegans/results/20191217/markers_to_5cm_bins.RDS')
cisMarkers=readRDS('/data/single_cell_eQTL/elegans/results/20191217/closestMarker.RDS')
G=standardise2(Gr)
#---------------------------------------------------------------------------------------------
Gcor=crossprod(G)/(nrow(G)-1)
fcm=caret::findCorrelation(Gcor,cutoff=.99999, verbose=F)
m.to.keep=which(!colnames(Gr)%in%colnames(Gr)[fcm])
Greduced=Gr[,m.to.keep]

lGreduced=log(Greduced/(1-Greduced))
lGreduced.var=colVars(lGreduced)
lGreduced=lGreduced[,which(!is.na(lGreduced.var))]
Greduced=Greduced[,which(!is.na(lGreduced.var))]

rm(Gcor)
m.to.keep=match(colnames(Greduced), colnames(Gr))

markerGRr=markerGR[m.to.keep]
tag=transcript.annotations[transcript.annotations$type=='gene',]
tag=tag[match(colnames(Yr), tag$gene_id),]
cisMarkers=data.frame(transcript=colnames(Yr), marker=colnames(Greduced)[nearest(tag,markerGRr)],stringsAsFactors=F)
cisMarkers=na.omit(cisMarkers)


# for each transcript, count up the number of cells per tissue the transcript is expressed in with at least one read
pte=apply(Yr, 2, function(x) sapply(split(x, joint.covariates$broad_tissue), function(y) sum(y>0)))
# hard cutoff for 'is expressed in' each broad tissue
t.cutoff=round(table(joint.covariates$broad_tissue)*.1)
t.cutoff20=20
is.expressed10=apply(pte,2,function(x) x>t.cutoff)
is.expressed20=apply(pte,2,function(x) x>t.cutoff20)
saveRDS(is.expressed10, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_10_percent_of_broad_tissue.RDS')
saveRDS(is.expressed20, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_20_cells_broad_tissue.RDS')

is.expressed20=readRDS('/data/single_cell_eQTL/elegans/results/20191217/expressed_in_20_cells_broad_tissue.RDS')


# str(cP)
# cPf=cP[cP$FDR<.05,]
# utt=unique(as.character(cPf$transcript))
# out of 6156 transcripts with sig linkages 271 fall in set of ~1000 transcripts with less than 20 cells with read counts in any tissue
# let's filter those 1000 transcripts!
expressed.transcripts=names(which(colSums(is.expressed20)>0))

# additional filters
Yr=Yr[,expressed.transcripts]
cisMarkers=cisMarkers[cisMarkers$transcript %in% expressed.transcripts,]

#tissueHotspotBetas=list()


#negative binomial null models within each tissue 
tissueNBNulls=list()


#functions to map LOD to FDR within each tissue 
qtl.fdrs=list()
# QTL peaks with FDR control within each tissue
tissueQTLPeaks=list()
# same as tissueQTLPeaks but with negative binomial regression betas, se, and p-vals for FDR<.2
tissueQTLPeaksNB=list()
# hotspots detected within each tissue
withinTissueHotspots=list()
# betas for hotspots detected within each tissue (all expressed transcripts in that tissue)
withinTissueHotspotBetas=list()
# betas for hotspots detected in the joint analysis 
jointTissueHotspotBetas=list()
# within tissue analysis for closest marker to each transcript
tissueCisBetas=list()

LODmatrix=list()
LPmatrix=list()
af.tissue=list()
tot.count.LOD=list()

fc=unique(as.character(joint.covariates$fine_tissue))
fc=fc[!is.na(fc)]
bc=unique(as.character(joint.covariates$broad_tissue))
bc=bc[!is.na(bc)]
cty=c(rep('broad', length(bc)), rep('fine', length(fc)))
types=c(bc,fc)
  
#20:86) { # 1:19 ) { #length(types)) {


# Workflow is run within tissue analysis then re-run hotspot fitting for joint analysis 
for(kk in 1:19) { 
    print(kk)
    uc=types[kk]
    print(uc)
    
    # tk is subset of cells classified as tissue =========================================
    if(cty[kk]=='broad') {
        tk=which(joint.covariates$broad_tissue==uc)
        dout=paste0('/data/single_cell_eQTL/elegans/results/20191217/broad/', uc, '/')
        print('broad')
    }
    if(cty[kk]=='fine') {
        tk=which(joint.covariates$fine_tissue==uc & !is.na(joint.covariates$PC1))
        dout=paste0('/data/single_cell_eQTL/elegans/results/20191217/fine/', uc, '/')
        print('fine')
    }
    dir.create(dout)
    print(kk)
    if(length(tk)<6) { next; }
    #======================================================================================

    
    #subset covariates ====================================================================
    covs=joint.covariates[tk,]
    #======================================================================================
   
   #plot(jitter(log2(colSums(Ns)),40),          jitter(log2(colSums(Cs)),40))

    # subset count matrix =================================================================
    Yk=Yr[tk,]
    tcounts=Yk>0
    tcounts=colSums(tcounts)
   
    # expressed in at least 20 cells  
    pick.filter=min(20, (length(tk)*.1))
    if(pick.filter<6) {pick.filter=6}
    Yk=Yk[,tcounts>pick.filter]
    #======================================================================================

    
    # subset genotype matrix ==============================================================
    Gsub=Greduced[tk,]#unique(cisMarkers$marker)]
  
    # base null model with batch and log2(total reads per cells) ==========================
    mmp1=model.matrix(lm(log(Yk[,1]+1)~log(covs$total)+covs$Batch))
    #======================================================================================

    #thought, fit cis-only model 
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% colnames(Yk)),]

    #svd on residuals
    # map that 
    nb.cisModel=fitNBCisModel(cMsubset, mmp1, Yk, Gsub)
    
    cisNB=nb.cisModel$cisNB
    nbR=nb.cisModel$nbR

    scresid1=svd(scale(nbR))
    plot(scresid1$d^2/sum(scresid1$d^2), xlim=c(0,20))
   
    PCs=scale((scresid1$u[,1:1000]))
    nt=as.integer(nrow(PCs))
    rff=crossprod(PCs,scale(Gsub))/(nt-1)
    LODr=-nt*log(1-rff^2)/(2*log(10))

    Yfe=scale(nbR) #)#Yresid2.list$nbR)
    Gf=scale(Gsub)
    #Gsubr=(residuals(lm(lGreduced[tk,]~log(covs$total)+covs$Batch)))
    #GsubrE=exp(Gsubr)/(1+exp(Gsubr))
    ##Gsubr.var=colVars(Gsubr)
    #Gf=scale(GsubrE)

    dataC=data.frame(cbind(mmp1, Gf))
    library(mpath)
    gn='WBGene00020612'
    test= glmreg(Yk[,gn]~., data=dataC, family="negbin", theta=cisNB[[gn]]$theta)

    nt=as.integer(nrow(Yfe))
    rff=crossprod(Yfe,Gf)/(nt-1)
    
    # don't retain standardized betas within tissue given difficulty 
    # of interpreting across tissues 
    #saveRDS(rff, paste0(dout, 'tissue_std_beta.RDS'))
    LODr=-nt*log(1-rff^2)/(2*log(10))
    LODmatrix[[uc]]=LODr
    
    ##saveRDS(fl(LODr), paste0(dout, 'tissue_LOD.RDS'))
    # and construct permutation nulls 
    nperm=20
    LODperm=replicate(nperm, {
                    nind=sample(1:nt)
                    pme=crossprod(Yfe[nind,],Gf)/(nt-1)
                    return((-nt*log(1-pme^2)/(2*log(10)))) 
                    })
    LPmatrix[[uc]]=LODperm

    tchr=transcript.data$chromosome_name[match(rownames(LODr),transcript.data$wbps_gene_id)]
    mchr=as.character(seqnames(mGRs))
    eg=expand.grid(tchr, mchr)
    be=which(eg$Var1==eg$Var2)

    #for(column in 1:ncol(emat)){emat[rownames(emat)==colnames(emat)[i],column]=1; print(column); }
    LODr[be]=0
    for(i in 1:nperm) {     LODperm[,,i][be]=0    }


    m1=which(mchr=='I')
    pmean=apply(LODperm[,m1,],c(2,3), mean)
    plot(colMeans(LODr[,m1]))

    mGRs=markerGR[match(colnames(Gf), colnames(Gr))]
    twgFDR= getFDRfx(LODr,LODperm, cisMarkers) 
    tcP=getGlobalPeaks(rff, LODr,twgFDR, mGRs,  transcript.data, chroms) 

    tcPsig=fitNBTransModel(tcP, Yk,Gsub, cisNB, cMsubset, .05) 
    tcPsig$FDR.negbin=p.adjust(tcPsig$negbin.p, 'fdr')
    t.hotspots = getHotspots(tcPsig, gmapd, fdr.thresh=.05, pFDR=T)
    ggplot(data.frame(rbindlist(t.hotspots[[2]])), aes(x=values, y=lengths, color=peaks))+geom_col()+
        facet_grid(~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('10cm bins across genome')

    qm=tcPsig[tcPsig$FDR<.01, ] #& tcPsig$FDR.negbin<.05,]
    ggplot(qm,aes(x=pos,y=tpos,alpha=-log10(FDR.negbin)))+geom_point()+
        geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) + #, color ="grey70", alpha=.1))+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') +  
        #geom_vline(aes(xintercept=start),data=wthl2,color='red',alpha=.75)+
        facet_grid(tchr~chr,scales="free")+theme_classic()  













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


   
#save results 
perTissueResults=list(
    qtl.fdrs=qtl.fdrs,
    tissueQTLPeaks=tissueQTLPeaks,
    tissueQTLPeaksNB=tissueQTLPeaksNB,
    withinTissueHotspots=withinTissueHotspots,
    withinTissueHotspotBetas=withinTissueHotspotBetas,
    jointTissueHotspotBetas=jointTissueHotspotBetas,
    tissueCisBetas=tissueCisBetas,
    tissueNBNulls=tissueNBNulls)
#saveRDS(perTissueResults, file='/data/single_cell_eQTL/elegans/results/20191217/perTissueLists.RDS')
perTissueResults=readRDS('/data/single_cell_eQTL/elegans/results/20191217/perTissueLists.RDS')


#==== Joint analysis ===============================================================================
# joint LOD is sum LOD scores for each transcript
# null comes from sum of permutations across tissues 
# pooling transcripts based on number of cells a transcript is classified as expressed in 

# get set of all transcripts with mapping stats 
utranscripts=unique(do.call('c', sapply(LODmatrix, function(x) rownames(x))))

# sum lods 
LODsums=matrix(0,length(utranscripts), ncol(LODmatrix[[1]])) 
rownames(LODsums)=utranscripts
colnames(LODsums)=colnames(LODmatrix[[1]])
for(tt in names(LODmatrix)){
    print(tt)
       x=LODmatrix[[tt]]
       LODsums[rownames(x),]=LODsums[rownames(x),]+x
}

# sum perm lods 
LODpsums=array(0,dim=c(length(utranscripts), ncol(LODmatrix[[1]]),10) )
rownames(LODpsums)=utranscripts
colnames(LODpsums)=colnames(LODmatrix[[1]])
for(tt in names(LPmatrix)){
    print(tt)
       x=LPmatrix[[tt]]
       LODpsums[rownames(x),,]=LODpsums[rownames(x),,]+x
}
#saveRDS(LODsums, file='/data/single_cell_eQTL/elegans/results/20191217/jLODsums.RDS')
#saveRDS(LODpsums, file='/data/single_cell_eQTL/elegans/results/20191217/jLODpermsums.RDS')

# total number of cells a transcript is expresssed in 
csYr=colSums( as.numeric(table(joint.covariates$broad_tissue))*is.expressed20)
csYr=csYr[match(colnames(Yr), colnames(is.expressed20))]
csYr=csYr[which(names(csYr)%in%utranscripts)]
#cis20s=split(colnames(is.expressed20), colSums(is.expressed20))
#csYr=colSums(Yr>0)
hist(csYr, breaks=100)
abline(v=Hmisc::cut2((csYr),m=1000, onlycuts=T))
# bin into ~ bins of ~1000 transcripts 

csYrc=as.character(Hmisc::cut2((csYr),m=1000))
cis20s=split(names(csYr),csYrc)

# for each bin get separate FDR function =================================================
jointPeaks=list()
jointCis=list()
for(cg in names(cis20s)) {
#cg="19"
    print(cg)
    L0=LODsums[cis20s[[cg]],]
    LPP=LODpsums[cis20s[[cg]],,]
    LG1f=getFDRfx(L0, LPP, cisMarkers)

    mGRs=markerGR[match(unique(cisMarkers$marker), colnames(Gr))]
    jointPeaks[[cg]]=getGlobalPeaks(L0, L0,LG1f, mGRs, transcript.data, chroms) 

    cMsubsetJ=cisMarkers[which(cisMarkers$transcript %in% rownames(L0)),]
    cisBetasJ=L0[cbind(cMsubsetJ$transcript, cMsubsetJ$marker)]
    names(cisBetasJ)=cMsubsetJ$transcript
    max.obsLODc=L0[cbind(cMsubsetJ$transcript, cMsubsetJ$marker)]

    cisBetasJ=data.frame(transcript=names(cisBetasJ), 
                           cisMarker=cMsubsetJ$marker,
                           #sbeta=cisBetasJ,
                           LOD=max.obsLODc,
                           FDR=LG1f[[4]](max.obsLODc), stringsAsFactors=F)
    cisBetasJ$FDR[is.na(cisBetasJ$FDR)]=1
    cisBetasJ$FDR[cisBetasJ$FDR>1]=1
    jointCis[[cg]]=cisBetasJ

}
#=============================================================================================

#jPs=lapply(jointPeaks, function(x) x[x$FDR<.1,])

#sapply(jPs, function(x) min(x$LOD))
jPs=rbindlist(jointPeaks,idcol='ntissues')
saveRDS(jPs, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_peaks.RDS')

jointHotspots=getHotspots(jPs, gmapd, fdr.thresh=.05, pFDR=T)
saveRDS(jointHotspots, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_hotspots.RDS')

#sapply(jPs, function(x) min(x$LOD))
cPs=rbindlist(jointCis,idcol='ntissues')
saveRDS(cPs, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_cis.RDS')

library(pheatmap)

#visualizations 
cPf=jPs[jPs$FDR<.01,]
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



for( tt in names(qm2)){
    qm=data.frame(qm2[[tt]])
    qm=qm2
    #wthl2=wthl[[tt]]
    #wthl2$chr=factor(wthl2$seqnames,levels=rev(c("X","V","IV","III","II","I")))
   
   
    ggplot(qm,aes(x=pos,y=tpos,alpha=-log10(FDR+1e-6)/6))+geom_point()+
        geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) + #, color ="grey70", alpha=.1))+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') +  
        #geom_vline(aes(xintercept=start),data=wthl2,color='red',alpha=.75)+
        facet_grid(tchr~chr,scales="free")+theme_classic()+
        theme_classic()+
        ggtitle(tt)
    ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/', tt,'_glm_nh.png'),width=10,height=10)
}
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



