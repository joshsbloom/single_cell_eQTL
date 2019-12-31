library(data.table)
library(ggsci)
library(ggplot2)
library(pracma)


#fdrfx.fc=getFDRfx(rff,Yfe,Gf,nperm=5,vint=seq(0.001,.45,.001))

generatePlot = function(jointPeaks,dfReturn=F,cutoff=0.1)
{
  colnames(jointPeaks)[colnames(jointPeaks)=="peak.marker"] = "marker"
  jointPeaks = data.frame(jointPeaks,chr=gsub("_+.*","",jointPeaks$marker),pos=as.numeric(gsub(".*_","",jointPeaks$marker)))
  jointPeaksGeneAnnots = geneAnnotsCur[match(jointPeaks$transcript,geneAnnotsCur$wbps_gene_id),]
  jointPeaks = jointPeaks[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  jointPeaksGeneAnnots = jointPeaksGeneAnnots[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  plot_df = data.frame(jointPeaks,jointPeaksGeneAnnots) %>% filter(FDR<cutoff)
  plot_df$chromosome_name = factor(plot_df$chromosome_name,levels=c("X","V","IV","III","II","I"))
  if(dfReturn)
  {
    plot_df
  }else
  {
 (ggplot(plot_df,aes(x=pos,y=start_position)) + facet_grid(chromosome_name~chr,scales="free") + geom_point(aes(alpha=-log10(FDR))) + theme_classic())
  }
}


getFDRfx=function(rff, Yfe,Gf,cisMarkers, nperm=5){

    max.obsLOD=rowMaxs(abs(rff), value=T) #apply(abs(rff),1,max)
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(rff)),]
    
    max.obsLODc=abs(rff[cbind(cMsubset$transcript, cMsubset$marker)])
    vint=seq(0.001,max(max.obsLOD)+.001, .001)
    
    # no point keeping all this, just retain max stats 
    ll=replicate(nperm, {
                    nind=sample(1:nrow(Yfe))
                    pme=abs(crossprod(Yfe[nind,],Gf)/(nrow(Yfe)-1))
                    lout=list(amax=rowMaxs(pme, value=T),
                              cmax=pme[cbind(cMsubset$transcript, cMsubset$marker)] 
                    )
                    return(lout)
        }, simplify=F )

  
    # global FDR  ---------------------------------------------------------------
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = vint

    #ll2=ll[c(2,4)]
    # expected number of QTL peaks with LOD greater than threshold
    ll1=lapply(ll, function(x) x$amax)
    #expPcnt = sapply(vint,  
    #             function(thresh) { 
    #                 return(mean(apply(ll1,2,function(x) sum(x>thresh))))
    #             })
    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(sapply(ll1,function(x) sum(x>thresh))))
    })

    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    
    pFDR = rev(cummax(rev(pFDR)))
    pFDR[1]=pFDR[1]-(1e-3)
    g.fdrFX=approxfun(pFDR, vint, ties=mean)
    g.rtoFDRfx=approxfun(vint,pFDR, ties=mean)


    # local FDR --------------------------------------------------------------------
    ll2=lapply(ll, function(x) x$cmax)
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLODc>thresh) }   )
    names(obsPcnt) = vint
    
    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(sapply(ll2,function(x) sum(x>thresh))))
                 })
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    
    pFDR = rev(cummax(rev(pFDR)))
    pFDR[1]=pFDR[1]-(1e-3)
    c.fdrFX=approxfun(pFDR, vint, ties=mean)
    c.rtoFDRfx=approxfun(vint,pFDR, ties=mean)
    #---------------------------------------------------------------------------------

    return(list(g.fdrFX=g.fdrFX,g.rtoFDRfx=g.rtoFDRfx,
                c.fdrFX=c.fdrFX, c.rtoFDRfx=c.rtoFDRfx))
}



# get peaks from genomewide scan ------------------------------------------------------------
getGlobalPeaks=function(rff, LODr,fdrfx.fc,markerGR,transcript.data, chroms) {
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
                             FDR=fdrfx.fc[['g.rtoFDRfx']](abs(mstat)) )
        cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
        cPeaksT[[cc]]$FDR[(cPeaksT[[cc]]$FDR)>1]=1
    }
 
    cP=rbindlist(cPeaksT, idcol='chrom')
    cP$peakGpos=markerGR$gcoord[match(cP$peak.marker, colnames(rff))]
    cP$tGpos=transcript.data$gcoord[match(cP$transcript, transcript.data$wormbase_gene)]
    cP$peakLGpos=markerGR$gcoord[match(cP$CI.l, colnames(rff))]
    cP$peakRGpos=markerGR$gcoord[match(cP$CI.r, colnames(rff))]
    #saveRDS(cP, file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_jointPeaks_filt.RDS')
    #cP=readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_jointPeaks_filt.RDS')
    cP$tchr=transcript.data$chromosome_name[match(cP$transcript, transcript.data$wormbase_gene)]
    cP$chr=tstrsplit(cP$peak.marker, '_')[[1]]
    cP$pos=as.numeric(tstrsplit(cP$peak.marker, '_')[[2]])
    return(cP)
}
#---------------------------------------------------------------------------------------------------





getHotspots=function(tcP, gmapd, fdr.thresh=.05) {

    ####### identify hotspots 
    cPsig=tcP[tcP$FDR<fdr.thresh,]
    cPsig.nocis=cPsig[cPsig$tchr!=cPsig$chr,]
    cPsig.nocis$cm5bin=gmapd$cm5bin[match(cPsig.nocis$peak.marker, gmapd$marker)]
    
    bcount=data.frame(values=1:130, lengths=0)

    bcountr=rle(sort(cPsig.nocis$cm5bin))
    bcount$lengths[bcountr$values]=bcountr$lengths
    #bcount=data.frame(values=bcount$values, lengths=bcount$lengths)
    #bcount.holder$lengths[

    sig.hp=qpois(1-(.05/max(gmapd$cm5bin)),ceiling(mean(bcount$lengths)))+1 
    gcm=t(sapply(gmapdc, function(x) range(x$cm5bin)))
    bcount$chr=rep(chroms, (1+gcm[,2]-gcm[,1]))
    bcount$sig=bcount$lengths>sig.hp

    bcounts=split(bcount,bcount$chr)

    #automate peak detection
    bcl=lapply(bcounts, function(x) {
           f=findpeaks(x$lengths,minpeakheight=sig.hp, minpeakdistance=3, nups=0)[,2]
           x$peaks=FALSE
           x$peaks[f]=T
           return(x)  })
    hotspot.bins=unlist(sapply(bcl, function(x) x$values[x$peaks]))
    cpb=cPsig.nocis[cPsig.nocis$cm5bin %in% as.vector(unlist(hotspot.bins)), ]
    cpbs=split(cpb, cpb$cm5bin)
    hotspot.markers=as.character(sapply(cpbs, function(x) names(which.max(table(x$peak.marker)))))
    return(list(hotspot.markers=hotspot.markers, bin.counts=bcl))
}

mvn.scanone=function(g.s, Y,add.cov=NULL, roi1=NULL) {
           ldetRSSf.1=rep(NA, ncol(g.s))
           if(!is.null(roi1)) { ss=roi1 } else { ss=1:ncol(g.s) } 
           for(i in min(ss):max(ss)) {
               #ldetRSSf.1[i]=det_AtA(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               if(is.null(add.cov) ) {
                RSSf.1=crossprod(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               } else {
                RSSf.1=crossprod(residuals(.lm.fit(cbind(add.cov, g.s[,i]), Y)))
               }
               ldetRSSf.1[i]=determinant(RSSf.1, logarithm=T)$modulus
           }
           return(ldetRSSf.1)
}


# cells with NAs are useless
#crit=!is.na(fine.classification$broad) #& !is.na(fine.classification$fine) 
# & log2(fine.classification$total)>8.5 & log2(fine.classification$total)<15.5

crit=!is.na(joint.covariates$broad_tissue)
joint.covariates=joint.covariates[crit,]
saveRDS(joint.covariates, '/data/single_cell_eQTL/elegans/results/20191217/covariates.RDS')

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

#geno=do.call('cbind', geno)
Gr=geno[crit,]
saveRDS(Gr, '/data/single_cell_eQTL/elegans/results/20191217/geno_matrix.RDS')

G=standardise2(Gr)
grs=split(data.frame(Gr), joint.covariates$broad_tissue)
gvar=lapply(grs, function(x) colVars(data.matrix(x), parallel=T))
saveRDS(gvar, '/data/single_cell_eQTL/elegans/results/20191217/broad/MarkerVariance.RDS')

#correlate genotypes across strains 
#A=tcrossprod(G)
#A=A/(nrow(G)-1)
#rr=which(A>.99, arr.ind=T)
#rr=rr[rr[,1]!=rr[,2],]

colnames(G)=colnames(G)
rownames(G)=rownames(Yr)

# bin the genome into 5cm bins -----------------------------------------------------------------------------
gmapd=gmap[match(colnames(G), gmap$marker),]
gmapdc=split(gmapd, gmapd$chrom)
gmap.max=sapply(gmapdc, function(x) ceiling(max(x$map)))
mbin=0
bin.size=5
for(cc in chroms) {
    fib=findInterval(gmapdc[[cc]]$map,seq(0,gmap.max[[cc]]+(bin.size-gmap.max[[cc]]%%bin.size),bin.size))
    gmapdc[[cc]]$cm5bin=fib+mbin
    mbin=max( gmapdc[[cc]]$cm5bin )
}
gmapd=rbindlist(gmapdc)
saveRDS(gmapd, file = '/data/single_cell_eQTL/elegans/results/20191217/markers_to_5cm_bins.RDS')
#---------------------------------------------------------------------------------------------------------

m0=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch)
r0=(residuals(lm.fit(m0,log2(Yr+1))))
m1=model.matrix(Yr[,1]~joint.covariates$broad_tissue-1)
ctm=lm.fit(m1, r0)
ctm=coef(ctm)
rownames(ctm)=gsub('joint.covariates\\$broad_tissue', '', rownames(ctm))
saveRDS(ctm, file = '/data/single_cell_eQTL/elegans/results/20191217/transcript_tissue_betas.RDS')

#model total reads, batch, broad and fine tissue classifications 
mm=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue) 
#+joint.covariates$fine_tissue)
BOLS=lm.fit(mm, log2(Yr+1), intercept=F)
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

Yresid=residuals(BOLS)
rownames(Yresid)=rownames(Yr)
colnames(Yresid)=colnames(Yr)
rm(BOLS)
#Yresid.var=colVars(Yresid,parallel=T)

Yre=standardise2(Yresid)
rownames(Yre)=rownames(Yresid)
colnames(Yre)=colnames(Yresid)

# identify closest marker to transcript for cis only test (to do) --------------------------------------------
g2m=tstrsplit(colnames(G),'_')
markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*",name=colnames(G)))
markerGR$gcoord=gcoord.key[as.character(seqnames(markerGR))]+start(markerGR)

tag=transcript.annotations[transcript.annotations$type=='gene',]
tag=tag[match(colnames(Yre), tag$gene_id),]
cisMarkers=data.frame(transcript=colnames(Yre), marker=colnames(G)[nearest(tag,markerGR)],stringsAsFactors=F)
cisMarkers=na.omit(cisMarkers)
saveRDS(cisMarkers, '/data/single_cell_eQTL/elegans/results/20191217/closestMarker.RDS')


#indices of cis eQTL in 
#rows are transcripts, columns are markers
cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
#------------------------------------------------------------------------------------------------------------

#sinfo=seqinfo(markerGR)
#sinfo@seqlengths=sapply(sapply(sapply(split(markerGR, seqnames(markerGR)),function(x) start(x)), max), as.integer)


# 1D scan speedup 
r=crossprod(Yre,G)/(nrow(Yre)-1)
saveRDS(r, '/data/single_cell_eQTL/elegans/results/20191217/joint/joint_std_beta.RDS')

LOD=-nrow(Yre)*log(1-r^2)/(2*log(10))
#saveRDS(r, '/data/single_cell_eQTL/elegans/results/joint_all_betas.RDS')
#saveRDS(LOD, '/data/single_cell_eQTL/elegans/results/joint_all_LODs.RDS')

#pool for permutations
wgFDR= getFDRfx(r, Yre, G, cisMarkers, nperm=5)
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

# inputs 1) all peaks, 2) markers sorted into centimorgan bins, 3, fdr threshold
hotspots = getHotspots(cP, gmapd)
saveRDS(hotspots,'/data/single_cell_eQTL/elegans/results/20191217/joint/hotspots.RDS')

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

#error on 45
for(kk in 20:86) { # 1:19 ) { #length(types)) {
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

    # inputs 1) all peaks, 2) markers sorted into centimorgan bins, 3, fdr threshold
    #tryCatch({
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



