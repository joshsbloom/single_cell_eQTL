install.packages('pracma')
install.packages('Rfast2')
install.packages('fastglm')
install.packages('float')

library(data.table)
library(ggsci)
library(ggplot2)
library(pracma)
library(Rfast2)
library(Rfast)
library(fastglm)
library(float)
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
gmapd=binMarkersToGeneticMap(Gr,gmap)
saveRDS(gmapd, file = '/data/single_cell_eQTL/elegans/results/20191217/markers_to_5cm_bins.RDS')
#-------------------------------------------------------------------------------------------------------

# identify closest marker to transcript for cis only test (to do) ----------------------------------
g2m=tstrsplit(colnames(G),'_')
# would functionalize this, but using markerGR later
markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*",name=colnames(G)))
markerGR$gcoord=gcoord.key[as.character(seqnames(markerGR))]+start(markerGR)

tag=transcript.annotations[transcript.annotations$type=='gene',]
tag=tag[match(colnames(Yr), tag$gene_id),]
cisMarkers=data.frame(transcript=colnames(Yr), marker=colnames(G)[nearest(tag,markerGR)],stringsAsFactors=F)
cisMarkers=na.omit(cisMarkers)
saveRDS(cisMarkers, '/data/single_cell_eQTL/elegans/results/20191217/closestMarker.RDS')
#-------------------------------------------------------------------------------------------------

#per tissue betas  (no QTL) -----------------------------------------------------------------------
m0=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch)
r0=(residuals(lm.fit(m0,log2(Yr+1))))
m1=model.matrix(Yr[,1]~joint.covariates$broad_tissue-1)
ctm=lm.fit(m1, r0)
ctm=coef(ctm)
rownames(ctm)=gsub('joint.covariates\\$broad_tissue', '', rownames(ctm))
saveRDS(ctm, file = '/data/single_cell_eQTL/elegans/results/20191217/transcript_tissue_betas.RDS')
#---------------------------------------------------------------------------------------------------


#if interrupted before line 82-----------------------------------------------------------------
joint.covariates=readRDS( '/data/single_cell_eQTL/elegans/results/20191217/covariates.RDS')
Yr=readRDS('/data/single_cell_eQTL/elegans/results/20191217/transcript_matrix.RDS')
Gr=readRDS('/data/single_cell_eQTL/elegans/results/20191217/geno_matrix.RDS')
gmapd=readRDS('/data/single_cell_eQTL/elegans/results/20191217/markers_to_5cm_bins.RDS')
cisMarkers=readRDS('/data/single_cell_eQTL/elegans/results/20191217/closestMarker.RDS')
G=standardise2(Gr)
#---------------------------------------------------------------------------------------------


# for each transcript, count up the number of cells per tissue the transcript is expressed in with at least one read
pte=apply(Yr, 2, function(x) sapply(split(x, joint.covariates$broad_tissue), function(y) sum(y>0)))
# hard cutoff for 'is expressed in' each broad tissue
t.cutoff=round(table(joint.covariates$broad_tissue)*.1)
t.cutoff20=20
is.expressed10=apply(pte,2,function(x) x>t.cutoff)
is.expressed20=apply(pte,2,function(x) x>t.cutoff20)
saveRDS(is.expressed10, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_10_percent_of_broad_tissue.RDS')
saveRDS(is.expressed20, file='/data/single_cell_eQTL/elegans/results/20191217/expressed_in_20_cells_broad_tissue.RDS')

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


fc=unique(as.character(joint.covariates$fine_tissue))
fc=fc[!is.na(fc)]
bc=unique(as.character(joint.covariates$broad_tissue))
bc=bc[!is.na(bc)]
cty=c(rep('broad', length(bc)), rep('fine', length(fc)))
types=c(bc,fc)
    #20:86) { # 1:19 ) { #length(types)) {

#LODmatrix=list()
#LPmatrix=list()

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
    Gsub=Gr[tk,unique(cisMarkers$marker)]
    #af.tissue[[uc]]=colSums(Gsub)/nrow(Gsub)
    #======================================================================================


    # base null model with batch and log2(total reads per cells) ==========================
    mmp1=model.matrix(lm(log2(Yk[,1]+1)~log2(covs$total)+covs$Batch))
    #======================================================================================


    #run once =============================================================================
    # goal here is to estimate dispersion within each tissue for each transcript
    #nbM0=negbinNull(colnames(Yk),Yk,mmp1)
    #saveRDS(nbM0, paste0(dout,'nb_null.RDS'))
    #tissueNBNulls[[uc]]=nbM0
    nbM0=tissueNBNulls[[uc]]
    #======================================================================================
    
    ## recalculate betas at joint hotspots, requires defining jointHotspots[[1]] ==========
    ## which happens outside this loop
    #jnbHotspots=negbinRefitHotspots(jointHotspots[[1]], Yk, Gsub,mmp1,nbM0)
    #saveRDS(jnbHotspots,paste0(dout, 'tissue_joint_hotspotBetas.RDS'))
    #jointTissueHotspotBetas[[uc]]=jnbHotspots
    ##=====================================================================================

   
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


    # for OLS Log2(Y+1) model get standardized genotypes ==================================    
    Gf=standardise2(Gsub)
    colnames(Gf)=colnames(Gsub)
    rownames(Gf)=rownames(Yk)
    #======================================================================================

    # Fast matrix-based mapping ===========================================================
    nt=as.integer(nrow(Yfe))
    # then map within each classification
    rff=crossprod(Yfe,Gf)/(nt-1)
    # don't retain standardized betas within tissue given difficulty 
    # of interpreting across tissues 
    #saveRDS(rff, paste0(dout, 'tissue_std_beta.RDS'))
    LODr=-nrow(Yfe)*log(1-rff^2)/(2*log(10))
    #LODmatrix[[uc]]=LODr
    ##saveRDS(fl(LODr), paste0(dout, 'tissue_LOD.RDS'))
    # and construct permutation nulls 
    #nperm=10
    #LODperm=replicate(nperm, {
    #                nind=sample(1:nrow(Yfe))
    #                pme=crossprod(Yfe[nind,],Gf)/(nt-1)
    #                return((-nrow(Yfe)*log(1-pme^2)/(2*log(10)))) 
    #                })
    #LPmatrix[[uc]]=LODperm
    ### faster to keep in memory than save and reload 
    ###saveRDS(apply(LODperm,3,fl), paste0(dout, 'tissue_LOD_perm.RDS'))
    #======================================================================================


    # return functions to map LOD to FDR for global and cis-only tests=====================
    #twgFDR= getFDRfx(LODr,LODperm, cisMarkers) 
    #qtl.fdrs[[uc]]=twgFDR
    #saveRDS(twgFDR, paste0(dout, 'FDRfx.RDS'))
    twgFDR=qtl.fdrs[[uc]]
    #======================================================================================
    

    # get peaks from global analysis within each tissue ===================================
    # mGRs=markerGR[match(unique(cisMarkers$marker), colnames(Gr))]
    # tcP=getGlobalPeaks(rff, LODr,twgFDR, mGRs,
    #                    transcript.data, chroms) 
    # tissueQTLPeaks[[uc]]=tcP
    # saveRDS(tcP, paste0(dout,'Peaks.RDS')) 
    tcP=tissueQTLPeaks[[uc]]  
    #======================================================================================


    # fit negative binomial regression model for each global peak if FDR < .2 =============
    # tcP_nb=negbinRefit(tcP, Yk,Gsub, mmp1, nbM0, fdr.thresh=.2)
    # trcp1=rbindlist(tcP_nb,fill=T)
    # trcp1f=trcp1[!is.na(trcp1$negbin.p),]
    # trcp1f$FDR.negbin=p.adjust(trcp1f$negbin.p, method='fdr')
    # saveRDS(trcp1f, paste0(dout,'PeaksNB.RDS')) 
    # tissueQTLPeaksNB[[uc]]=trcp1f #tcP
    trcp1f= tissueQTLPeaksNB[[uc]] 
    #======================================================================================

    # analysis of hotspots within each tissue =============================================
    #t.hotspots = getHotspots(trcp1f, gmapd, fdr.thresh=.05, pFDR=F)
    t.hotspots= withinTissueHotspots[[uc]]
    # withinTissueHotspots[[uc]]=t.hotspots
    # saveRDS(t.hotspots,paste0(dout, 'within_tissue_hotspots.RDS'))
    
    # calc negative binomial model and within-tissue hotspots 
    if( length(t.hotspots[[1]])>0 ) {
        wtHotspots=negbinRefitHotspots(t.hotspots[[1]], Yk, Gsub,mmp1,nbM0)
        #withinTissueHB=rff[,tissue.hotspots[[1]]]
        #saveRDS(wtHotspots,paste0(dout, 'within_tissue_hotspot_betas.RDS'))
        withinTissueHotspotBetas[[uc]]=wtHotspots
    }
    #=======================================================================================


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
    cisNB=cisNegbinFit(cMsubset,Yk,Gsub,mmp1,nbM0)
    cisNB=rbindlist(cisNB,fill=T,idcol='transcript')
    cisBetasTis=merge(cisBetasTis, cisNB, by='transcript', all.x=T,sort=F)
    cisBetasTis$FDR.negbin=p.adjust(cisBetasTis$negbin.p, method='fdr')
    tissueCisBetas[[uc]]=cisBetasTis
    #saveRDS(cisBetasTis, paste0(dout, 'tissue_cisBetas.RDS'))
    #========================================================================================
}

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


#visualizations 
cPf=jPs[jPs$FDR<.01,]
cPf$alpha=(-log10(cPf$FDR+1e-6)/6)
cPf$gmap.pos=gmapd$map[match(cPf$peak.marker, gmapd$marker)]

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
cPf$alpha=-log10(cPf$FDR.negbin+1e-6)/6



