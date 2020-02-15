##BiocManager::install("qvalue")
##install.packages('float')
library(Matrix)
library(rtracklayer)
library(data.table)
library(parallel)
library(Rfast)
library(Rfast2)
library(data.table)

library(ggplot2)
library(ggsci)

library(pracma)

library(MASS)
library(mgcv)
library(fastglm)

library(caret)
library(parallel)
library(abind)

source('/data/single_cell_eQTL/elegans/code/mapping_fxs.R')

reference.dir='/data/single_cell_eQTL/elegans/reference/'
data.structures.dir='/data/single_cell_eQTL/elegans/results/20200201/'


#correlate genotypes across strains 
#A=tcrossprod(G)
#A=A/(nrow(G)-1)
#rr=which(A>.99, arr.ind=T)
#rr=rr[rr[,1]!=rr[,2],]
#---------------------------------------------------------------------------------------------------

#mapping.R can now be run separately from main.R
transcript.data=readRDS(paste0(reference.dir, 'transcriptData.RDS'))
transcript.annotations=import(paste0(reference.dir, 'genes.gtf'))
joint.covariates=readRDS(paste0(data.structures.dir, 'covariates.RDS'))
Yr=readRDS(paste0(data.structures.dir,'transcript_matrix.RDS'))
Gr=readRDS(paste0(data.structures.dir,'geno_matrix.RDS'))


# genetic map (this is for 10 generations) ------------------------------------
 gmap=readRDS(paste0(reference.dir, 'geneticMapXQTLsnplist.rds'))
 gmap$marker=paste0(gmap$chrom, '_', gmap$pos)
    #https://www.genetics.org/content/170/2/875
    #map size = X(AIL(J))=jX/2 (where j is number of generations); 
    #expected map size= 4*300/2 = 600cm; 1500cm/600cm = 2.5 (correction factor)
    # correction for F4
 gmap$map=gmap$map/2.5
#------------------------------------------------------------------------------

gmapd=binMarkersToGeneticMap(Gr,gmap, bin.size=30)

g2m=tstrsplit(colnames(Gr),'_')
# would functionalize this, but using markerGR later
markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*",name=colnames(Gr)))
#---------------------------------------------------------------------------------------------

#fast computation of marker correlations 
G=standardise2(Gr)
Gcor=crossprod(G)/(nrow(G)-1)
rm(G)
# visualize auto-correlation. e.g. correlation decay from given marker to a marker with an index 100 less
#v=Gcor[row(Gcor)==(col(Gcor)-100)] ... eek
# prune correlated markers  (approx every 5cm)
# previous code transformed genotypes as log(g/(1-g)) for a beta-regression-like setup
# ~30 or markers were lost by div by 0 errors, fixed 2/14
fcm=caret::findCorrelation(Gcor,cutoff=.9999, verbose=F)
m.to.keep=which(!colnames(Gr)%in%colnames(Gr)[fcm])
Greduced=Gr[,m.to.keep]
rm(Gcor)
m.to.keep=match(colnames(Greduced), colnames(Gr))

markerGRr=markerGR[m.to.keep]

tag=transcript.annotations[transcript.annotations$type=='gene',]
tag=tag[match(colnames(Yr), tag$gene_id),]
cisMarkers=data.frame(transcript=colnames(Yr), marker=colnames(Greduced)[nearest(tag,markerGRr)],stringsAsFactors=F)
cisMarkers=na.omit(cisMarkers)


## for each transcript, count up the number of cells per tissue the transcript is expressed in with at least one read
## hard cutoff for 'is expressed in' each broad tissue
pte=apply(Yr, 2, function(x) sapply(split(x, joint.covariates$broad_tissue), function(y) sum(y>0)))
t.cutoff20=20
is.expressed20=apply(pte,2,function(x) x>t.cutoff20)
#saveRDS(is.expressed20, file=paste0(data.strucutres.dir, 'expressed_in_20_cells_broad_tissue.RDS')) 
is.expressed20=readRDS(paste0(data.strucutres.dir, 'expressed_in_20_cells_broad_tissue.RDS'))

# would probably benefit from harsher filter here
expressed.transcripts=names(which(colSums(is.expressed20)>0))

# additional filters
Yr=Yr[,expressed.transcripts]
cisMarkers=cisMarkers[cisMarkers$transcript %in% expressed.transcripts,]

fc=unique(as.character(joint.covariates$fine_tissue))
fc=fc[!is.na(fc)]
bc=unique(as.character(joint.covariates$broad_tissue))
bc=bc[!is.na(bc)]
cty=c(rep('broad', length(bc)), rep('fine', length(fc)))
types=c(bc,fc)
  
# slow slow slow 
cl <- makeCluster(36)
clusterEvalQ(cl, {
       library(fastglm)
       library(Matrix)
       library(MASS) 
       NULL  })

#af.tissue=list()
#tot.count.LOD=list()

#functions to map LOD to FDR within each tissue 
qtl.fdrs=list()
# within tissue analysis for closest marker to each transcript
tissueCisNB=list()
# QTL peaks with FDR control within each tissue
tissueQTLPeaksNB=list()

# store results from LOD score calculation
LODmatrix=list()
# store results from permutation procedure
LPmatrix=list()

# within tissue analysis for closest marker to each transcript, separate structure for fine tissue categories
tissueCisNBFine=list()

dir.create(paste0(data.strucures.dir, '/broad/'))
dir.create(paste0(data.strucures.dir, '/fine/'))

# Map within each tissue  
for(kk in seq_along(cty)) { 
    print(kk)
    uc=types[kk]
    print(uc)
    

    # tk is subset of cells classified as tissue =========================================
    if(cty[kk]=='broad') {
        tk=which(joint.covariates$broad_tissue==uc)
        dout=paste0(data.structures.dir, '/broad/', uc, '/')
        print('broad')
    }
    if(cty[kk]=='fine') {
        tk=which(joint.covariates$fine_tissue==uc & !is.na(joint.covariates$PC1))
        dout=paste0(data.structures.dir, '/fine/', uc, '/')
        print('fine')
    }
    dir.create(dout)
    print(kk)
    if(length(tk)<6) { next; }
    #======================================================================================
    
    # tissueCisNB[[uc]]= readRDS(paste0(dout,'nb_cisModel.RDS'))
    # qtl.fdrs[[uc]]   = readRDS(paste0(dout, 'FDRfx.RDS'))
    # tissueQTLPeaksNB[[uc]]=readRDS(paste0(dout,'PeaksNB.RDS'))

    
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
    Gsub=Greduced[tk,]              #unique(cisMarkers$marker)]
  
    # base null model with batch and log2(total reads per cells) ==========================
    mmp1=model.matrix(lm(log(Yk[,1]+1)~log(covs$total)+covs$Batch))
    #======================================================================================

    #thought, fit cis-only model 
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% colnames(Yk)),]

  
    #negbin local marker test only
    cisNB=fitNBCisModel(cMsubset, mmp1, Yk, Gsub)          #, resids='deviance')
    if(cty[kk]=='broad') {
        saveRDS(cisNB, paste0(dout,'nb_cisModel.RDS'))
        tissueCisNB[[uc]]=cisNB
     }
    if(cty[kk]=='fine') {
        #if we're looking at a fine tissue then do not proceed with whole genome mapping
        tissueCisNBFine[[uc]]=cisNB
        next;
    }

    #cisNB=tissueCisNB[[uc]]

    thetas=sapply(cisNB, function(x) as.numeric(x$theta))
   
    Yks=Matrix(Yk, sparse=T)
    clusterExport(cl, varlist=c("thetas", "Yks", "Gsub", "mmp1", "nbLL", "domap"))
    clusterEvalQ(cl,{ Y=Yks;    DM=mmp1;   return(NULL);})
    
    #negbin QTL model
    LOD=do.call('rbind', parLapply(cl, names(thetas), domap) )
    LOD[is.na(LOD)]=0
    LOD=LOD/(2*log(10))
    rownames(LOD)=names(thetas)

    #permute within batch
    LODpL=list()
    for(i in 1:5) { 
        print(i)
        new.index=as.vector(do.call('c', sapply(split(seq_along(covs$Batch), as.character(covs$Batch)), sample)))
        Ykp=Matrix(Yk[new.index,], sparse=T)
        mmp2=mmp1
        mmp2[,2]=mmp1[new.index,2]
        clusterExport(cl, varlist=c('mmp2', 'Ykp'))
        clusterEvalQ(cl,{ Y=Ykp; DM=mmp2; return(NULL);})
        
        LODp=do.call('rbind', parLapply(cl, names(thetas), domap) )
        LODp[is.na(LODp)]=0
        LODp=LODp/(2*log(10))
        rownames(LODp)=names(thetas)
       
        LODpL[[as.character(i)]]=LODp
    }
    LODp=abind(LODpL, along=3)
    rm(LODpL)
    
    saveRDS(LOD, paste0(dout,'LOD.RDS')) 
    saveRDS(LODp, paste0(dout,'LODperm.RDS')) 

    LODmatrix[[uc]]=LOD
    LPmatrix[[uc]]=LODp

    twgFDR= getFDRfx(LOD,LODp, cMsubset) 
    saveRDS(twgFDR, paste0(dout, 'FDRfx.RDS'))
    qtl.fdrs[[uc]]=twgFDR
    #twgFDR=qtl.fdrs[[uc]]

    # sanity check
    #sum(twgFDR[[4]](sapply(cisNB,function(x)x$negbin.LRS/(2*log(10))))<.1)
    #sum(p.adjust(sapply(cisNB, function(x) x$negbin.p), method='fdr')<.1)
   
    mGRs=markerGR[match(colnames(Gsub), colnames(Gr))]
    tcP=getGlobalPeaks( LOD,twgFDR, mGRs,  transcript.data, chroms) 
    tcPsig=fitNBTransModel(tcP, Yk,Gsub, cisNB, mmp1, cMsubset, .5) 
    saveRDS(tcPsig, paste0(dout,'PeaksNB.RDS')) 
    tissueQTLPeaksNB[[uc]]=tcPsig 

}
stopCluster(cl)

perTissueList=list(
    qtl.fdrs=qtl.fdrs,
    # within tissue analysis for closest marker to each transcript
    tissueCisNB=tissueCisNB,
    # QTL peaks with FDR control within each tissue
    tissueQTLPeaksNB=tissueQTLPeaksNB)
saveRDS(perTissueList, paste0(data.structures.dir, "perTissueList.RDS"))


saveRDS(LODmatrix,paste0(data.structures.dir, "LODS.RDS")) 
saveRDS(LPmatrix, paste0(data.structures.dir, "LODS_5perms.RDS") )

#LODmatrix=readRDS("/data/single_cell_eQTL/elegans/results/20200201/LODS.RDS") 
#LPmatrix=readRDS("/data/single_cell_eQTL/elegans/results/20200201/LODS_5perms.RDS") 
#load in the per tissue lists
#attach(readRDS("/data/single_cell_eQTL/elegans/results/20200201/perTissueList.RDS"))


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
    LODpsums=array(0,dim=c(length(utranscripts), ncol(LODmatrix[[1]]),5) )
    rownames(LODpsums)=utranscripts
    colnames(LODpsums)=colnames(LODmatrix[[1]])
    for(tt in names(LPmatrix)){
        print(tt)
           x=LPmatrix[[tt]]
           LODpsums[rownames(x),,]=LODpsums[rownames(x),,]+x
    }

    # total number of cells a transcript is expresssed in 
    csYr=colSums( as.numeric(table(joint.covariates$broad_tissue))*is.expressed20)
    csYr=csYr[match(colnames(Yr), colnames(is.expressed20))]
    csYr=csYr[which(names(csYr)%in%utranscripts)]
    hist(csYr, breaks=100)
    abline(v=Hmisc::cut2((csYr),m=1000, onlycuts=T))
    # bin into ~ bins of ~1000 transcripts 

    csYrc=as.character(Hmisc::cut2((csYr),m=1000))
    cis20s=split(names(csYr),csYrc)

    # for each bin get separate FDR function =================================================
    jointPeaks=list()
    jointCis=list()
    for(cg in names(cis20s)) {
        print(cg)
        L0=LODsums[cis20s[[cg]],]
        LPP=LODpsums[cis20s[[cg]],,]
        LG1f=getFDRfx(L0, LPP, cisMarkers)

        mGRs=markerGR[match(unique(cisMarkers$marker), colnames(Gr))]
        jointPeaks[[cg]]=getGlobalPeaks(L0, LG1f, mGRs, transcript.data, chroms) 

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


# reorganize data structures and save --------------------------------------------------------
cPs=rbindlist(jointCis,idcol='ntissues')
saveRDS(cPs,paste0(data.structures.dir, 'jointModel_cisOnly.RDS'))

cisResults=lapply(tissueCisNB, 
                  function(x) data.frame(t(sapply(x, function(y)  do.call('c', y[1:8] ))),
                                         stringsAsFactors=F))
cisResults=lapply(cisResults, function(x) {
                      for( i in 3:8) {
                          x[,i]=as.numeric(x[,i]) }
                      return(x) })
for(tt in names(cisResults)){
    twgFDR=qtl.fdrs[[tt]]
    cisResults[[tt]]$FDR=twgFDR[[4]](cisResults[[tt]]$negbin.LRS/(2*log(10)))
    cisResults[[tt]]=cisResults[[tt]][!is.na(cisResults[[tt]]$FDR),]
}
saveRDS(cisResults, paste0(data.structures.dir,'perTissue_cisOnly.RDS')) 

jPs=rbindlist(jointPeaks,idcol='ntissues')
saveRDS(jPs, paste0(data.structures.dir,'jointModel_globalPeaks.RDS'))

saveRDS(tissueQTLPeaksNB, paste0(data.structures.dir,'perTissue_globalPeaks.RDS'))

jointHotspots=getHotspots(jPs, gmapd, fdr.thresh=.2)
saveRDS(jointHotspots, paste0(data.structures.dir,'jointModel_hotspots.RDS'))

wthl=lapply(tissueQTLPeaksNB, function(x) {
           getHotspots(x, gmapd, fdr.thresh=.2) })
saveRDS(wthl, paste0(data.structures.dir,'perTissue_hotspots.RDS'))




#Process results from fine tissue analysis ---------------------------------------------------
cisResultsFine=lapply(tissueCisNBFine, 
                  function(x) data.frame(t(sapply(x, function(y)  do.call('c', y[1:8] ))),
                                         stringsAsFactors=F))
cisResultsFine=lapply(cisResultsFine, function(x) {
                      for( i in 3:8) {
                          x[,i]=as.numeric(x[,i]) }
                      return(x) })
cisResultsFine=lapply(cisResultsFine, function(x) {
                      y=x
                      y$FDR=p.adjust(y$negbin.p, method='fdr')
                      y=na.omit(y)
                      return(y)
                      })

saveRDS(cisResultsFine, paste0(data.structures.dir. 'perFineTissue_cisOnly.RDS'))
#----------------------------------------------------------------------------------------------



# quick visualizations ------------------------------------------------------------------------
plotEQTL=function(cPf, titleme='') {
    ggplot(cPf,aes(x=pos,y=tpos,alpha=-log10(FDR+1e-10)/6))+geom_point()+
        geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) + 
        xlab('')+ylab('')+ scale_alpha(guide = 'none') +  
        facet_grid(tchr~chr,scales="free")+theme_classic()+ggtitle(titleme)
}

#of hotspots
wthsig=lapply(wthl, function(x) rbindlist(x[[2]]))
wthsig=wthsig[which(sapply(wthsig, function(x) sum(x$peaks))>0)]
wthsig$joint=rbindlist(jointHotspots[[2]])
hmerged=rbindlist(wthsig, idcol='tissue')

hmerged$tissue=relevel(as.factor(hmerged$tissue), ref="joint")
ggplot(hmerged, aes(x=values, y=lengths, color=peaks))+geom_col()+
    facet_grid(tissue~chr, scales="free")+ylab('trans-eQTL linkages')+xlab('30cm bins across genome')

#of joint eQTL

cPf=jPs[jPs$FDR<.2,]
plotEQTL(cPf, 'joint')

for( tt in names(tissueQTLPeaksNB)){
    qm=data.frame(tissueQTLPeaksNB[[tt]])
    #readline()
    plotEQTL(qm,tt)
    #ggsave(paste0('/home/jbloom/Dropbox/Lab Meeting - Presentations/2020/', tt,'_NB_glm_nh.png'),width=10,height=10)
}

