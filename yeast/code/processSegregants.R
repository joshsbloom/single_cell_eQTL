library(monocle)
library(Matrix)
library(rtracklayer)
library(data.table)
library(qtl)
library(qtl2)
library(parallel)
library(Rfast)
library(seqinr)
library(GenomicRanges)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(foreach)
library(doMC)
library(glmmTMB)
library(ggplot2)

# for parallelizing the HMM
ncores=16
registerDoMC(cores=ncores)

code.dir='/data/single_cell_eQTL/yeast/code/'
source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))


chroms=paste0('chr', as.roman(1:16)) 
#c('I','II', 'III', 'IV', 'V', 'X')

#some genome annotation
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

# for now load all crosses--------------------------------------------------
load(paste0(reference.dir, 'cross.list.RData'))
load(paste0(reference.dir, 'parents.list.RData'))
#reorganize
crossL=list('A'=cross.list[['A']], 
            'B'=cross.list[['B']], 
            '3004'=cross.list[['3004']],
            '393'=cross.list[['393']]

)
parentsL=list('A'=parents.list[['A']], 
              'B'=parents.list[['B']], 
              '3004'=parents.list[['3004']],
              '393'=parents.list[['393']]
)
rm(cross.list)
rm(parents.list)

# need genetic map
genetic.maps=lapply(crossL, extractGeneticMapDf)
#----------------------------------------------------------------------------

# set some variables---------------------------------------------------------
base.dir='/data/single_cell_eQTL/yeast/'
cranger.dir='filtered_feature_bc_matrix/'


cList=list(
           '01_2444_44_1-2'='B',
           '02_2444_44-'='B',
           '03_ByxRM_51-_1-2'='A',
           '04_ByxRM_51-'='A',
           '00_BYxRM_480MatA_1'='A',
           '00_BYxRM_480MatA_2'='A',
           '08_2444_cross_10k_Feb_21'='B',
           '09_2444_cross_5k_Feb_21'='B',
           '10_3004_cross_10k_Feb_21'='3004',
           '11_3004_cross_5k_Feb_21'='3004',
           '07_2444-cross-1'='B',
           '07_2444-cross-2'='B',
           '19_393_10k_May10'='393',
           '20_393_20k_May10'='393',
           '21_3004_10k_May10'='3004',
           '22_3004_20k_May10'='3004'
)

hList=list(
           '01_2444_44_1-2'=2^6.1,
           '02_2444_44-'= 2^7.5,
           '03_ByxRM_51-_1-2'=2^7,
           '04_ByxRM_51-'=2^7.4,
           '00_BYxRM_480MatA_1'=2^6.5,
           '00_BYxRM_480MatA_2'=2^6.5,
           '08_2444_cross_10k_Feb_21'=2^5.2,
           '09_2444_cross_5k_Feb_21'= 2^4.9,
           '10_3004_cross_10k_Feb_21'= 2^5.6,
           '11_3004_cross_5k_Feb_21'=  2^4.9,
           '07_2444-cross-1'  = 2^6,
           '07_2444-cross-2'  = 2^5.5,
           '19_393_10k_May10' = 2^5.2,
           '20_393_20k_May10' = 2^5.1,
           '21_3004_10k_May10'= 2^6,
           '22_3004_20k_May10'= 2^6.7
)


#experiments=c('1_2444_44_1-2', '2_2444_44-', '3_ByxRM_51-_1-2', 
#              '4_ByxRM_51-', '0_BYxRM_480MatA_1', '0_BYxRM_480MatA_2')
#cdesig=c('B','B', 'A', 'A', 'A', 'A')
#cell.cycle.classification.names=c('set1', 'set2', 'set3', 'set4', 'set01', 'set02')
#cell.cycle.classification.dir='results/cell_cycle/yeast_cc_annotations_v5/'


#experiments=c('08_2444_cross_10k_Feb_21',
#              '09_2444_cross_5k_Feb_21',
#              '10_3004_cross_10k_Feb_21',
#              '11_3004_cross_5k_Feb_21',
#              '07_2444-cross-1',
#              '07_2444-cross-2'
#        )
#cdesig=c('B', 'B',  '3004', '3004', 'B', 'B')            
              
              
experiments=names(cList)
cdesig=as.vector(sapply(cList, function(x) x))
het.thresholds=as.vector(sapply(hList, function(x) x))
data.dirs=paste0(base.dir, 'processed/', experiments, '/')
het.across.cells.threshold=.1
              
              
              
#              '1_2444_44_1-2', '2_2444_44-', '3_ByxRM_51-_1-2', 
#              '4_ByxRM_51-', '0_BYxRM_480MatA_1', '0_BYxRM_480MatA_2')
#cdesig=c('B','B', 'A', 'A', 'A', 'A')
#cell.cycle.classification.names=c('set1', 'set2', 'set3', 'set4', 'set01', 'set02')
#cell.cycle.classification.dir='results/cell_cycle/yeast_cc_annotations_v5/'
#by visual inspection (might want to revist this classification)
                # for 08-11,07
#het.thresholds=c(2^5.2, 2^4.9, 2^5.6, 2^4.9, 2^6,2^5.5)
                #for 1,2,3,4,0,0
                 #c(2^6.1, 2^7.5, 2^7, 2^7.4, 2^6.5,2^6.5)


#-----------------------------------------------------------------------------

sgd.genes=sgd.granges[sgd.granges$type=='gene',]
sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)

# iterate over each experiment -----------------------------------------------
for(ee in c(1:6,13:16) ) { #5:length(experiments)){
    experiment=experiments[ee]
    cross=crossL[[cdesig[ee]]]
    parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')
    dir.create(results.dir)

    #clean up this crap-------------------------
    gmap=genetic.maps[[cdesig[ee]]]
    gmap.s=gmap; colnames(gmap.s)[2]='map'
    gmap.s=split(gmap.s, gmap.s$chrom)
    gmap.s=gmap.s[chroms]
    gmap.s=jitterGmap(gmap.s)
    saveRDS(gmap.s, file=paste0(results.dir, 'gmap.RDS'))
    #-------------------------------------------
   
    #read in counts------------------------------------------------
        counts=readMM(paste0(data.dir,cranger.dir,'matrix.mtx.gz'))
        features=read.csv(paste0(data.dir,cranger.dir,'features.tsv.gz'), sep='\t',header=F, stringsAsFactors=F)
        barcodes=read.csv(paste0(data.dir,cranger.dir,'barcodes.tsv'),sep='\t',header=F,stringsAsFactors=F)[,1]
        # fix order of count matrix (put transcripts in genomic order) 
        reorderme=(order(match(features[,1], sgd$Name)))
        counts=counts[reorderme,]
        features=features[reorderme,]
        rownames(counts)=features[,1]
        colnames(counts)=barcodes
        saveRDS(counts, file=paste0(results.dir, 'counts.RDS'))
    #-------------------------------------------------------------

    # counts per parental allele (ref.counts and alt.counts, names not accurate) (fixes phasing) --------
    # also filters out counts from heavily distorted variants /// update this with diploid ASE data or parental data when we have it
        g.counts=getGenoInformativeCounts(data.dir, cranger.dir, gmap,parents, log2diff=8)
        saveRDS(g.counts, file=paste0(results.dir, 'gcounts.RDS'))
    #-----------------------------------------------------------------------------------------------------
    
    # can filter here for only sites that have informative variants in this data set,
    # but this may make comparison across datasets difficult in the future, hold on this for now
        rQTL.coded=encodeForRQTL(g.counts)
        #saveRDS(rQTL.coded, file=paste0(results.dir, 'rQTLcoded.RDS'))
    #----------------------------------------------------------------------------------------------------

    het.count=apply(rQTL.coded, 2, function(x) sum(x==0,na.rm=T))
    tot.count=colSums(g.counts[[1]]+g.counts[[2]]) #ref.counts+alt.counts)
   
    plot(log2(tot.count),log2(het.count),col="#00000022")  # manual inspection of total counts vs sites classified as hets
    
    # note this is classification as probably haploid (TRUE) or not (dipoid, remated segregant/doublet FALSE)
    classification=het.count<het.thresholds[ee]
    saveRDS(classification, file=paste0(results.dir, 'cellFilter.RDS'))

    png(file=paste0(results.dir, 'plot_classify_haploids.png'), width=1024, height=1024)
    ## add information about number of cells classified as haploid vs not
    plot(log2(tot.count),log2(het.count), col=  classification+1,
         xlab='log2(total geno informative reads)', ylab='log2(total heterozygous sites)', 
         main=experiment, sub=paste('total = ', length(het.count),  ' haploid =', sum(classification) ) )#s[[experiment]])))
    dev.off()

    # identify unreliable sites (sites called as hets in some fraction of cells, re-evaluate threshold chosen here)
    het.cells=!classification
    recurring.hets=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==0,na.rm=T)) > (het.across.cells.threshold*sum(!het.cells))
    print(sum(recurring.hets))
    print(sum(recurring.hets)/nrow(rQTL.coded)) 
    #per suggestion of James, run HMM on everything to start
    
    #cross2=buildCross2(rQTL.coded, gmap, 'riself', het.sites.remove=T, het.cells.remove=F,
    #                   het.cells=het.cells, recurring.hets=recurring.hets)
    #rqtl.genoprobs=calc_genoprob(cross2, error_prob=.005, lowmem=F, quiet=F, cores=16)
    #rqtl.genoprobs=do.call('cbind', sapply(rqtl.genoprobs, function(x) x[,2,]))
    #saveRDS(rqtl.genoprobs, file=paste0(results.dir, 'rqtl_genoprobs.RDS'))
    #emissionProbs=estimateEmissionProbs(g.counts[[1]],  g.counts[[2]], error.rate=.005,
    #                                    recurring.het.sites.remove=T,
    #                                    recurring.hets=recurring.hets)
    #test=runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms[3], calc='viterbi',return.chr=T) #, n.indiv=1000)
    #par(yaxs='i')
    #plot(colSums(test[which(classification),]-1)/sum(classification), ylim=c(0,1))
    #abline(h=.85)
    rm(rQTL.coded)
    emissionProbs=estimateEmissionProbs(g.counts[[1]],  g.counts[[2]], error.rate=.005,
                                        recurring.het.sites.remove=T,
                                        recurring.hets=recurring.hets)
    runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='genoprob') #, n.indiv=1000)
    runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='viterbi') #, n.indiv=1000)
}


#sgd.genes=sgd.granges[sgd.granges$type=='gene',]

# parallelized standardise function
standardise2 <- function (x, center = TRUE, scale = TRUE) {
  if ( center & scale ) {
    y <- t(x) - Rfast::colmeans(x,parallel=T)
    y <- y / sqrt( Rfast::rowsums(y^2) ) * sqrt( (dim(x)[1] - 1) )
	y <- t(y)
  } else if ( center & !scale ) {
    m <- Rfast::colmeans(x,parallel=T)
    y <- eachrow(x, m, oper ="-" )
  } else if ( !center & scale ) {
    s <- Rfast::colVars(x, std = TRUE,parallel=T)
    y <- eachrow(x, s, oper = "/")
  }
  colnames(y)=colnames(x)
  rownames(y)=rownames(x)
  y
} 

source(paste0(code.dir, 'ASE_fxs.R'))

library(mgcv)
library(Rfast2)
library(fastglm)
library(MASS)
library(irlba)
library(nebula)
library(dplyr)

fitNBCisModel=function(cMsubset, mmp1, Yk, Gsub, resids=NULL) {
    mmp0=mmp1[,-1]

    #nbR=Yk
    #nbR[is.numeric(nbR)]=NA
    ##nbP=nbR


    cisNB=list() 
    for(r in 1:nrow(cMsubset)){
        print(r)
        gn=as.character(cMsubset$transcript[r])
        p=as.character(cMsubset$marker[r])
        cmodel=cbind(mmp0, Gsub[,p])
        nbr=negbin.reg(Yk[,gn], cmodel,maxiters=500)
        if(!is.na(nbr$info[4]) & (nbr$info[4]<100 & !is.na(nbr$info[3])) ){
             theta.est=nbr$info[4]
             coef.est=nbr$be[,1]
        } else {
            nbr=mgcv::gam(Yk[,gn]~cmodel, family=mgcv::nb , method="ML")
            theta.est=nbr$family$getTheta(T)
            coef.est=as.vector(coef(nbr))
        }
        XX=cbind(mmp1, Gsub[,p])
        msize=ncol(XX)
        fnbrF=fastglm(XX, Yk[,gn], start=coef.est, family=negative.binomial(theta=theta.est,link='log'), maxit=500)
        ##pseudo r^2 statistics 
        ##library(performance)
        ##class(fnbrF)='glm'
        ##r2(fnbrF)
        fnbrN=fastglm(mmp1,Yk[,gn], start=coef.est[-msize],  family=negative.binomial(theta=theta.est,link='log'), maxit=500)
        #if(resids=='deviance'){       nbR[,gn]=residuals(fnbrF,'deviance') }
        #if(resids=='pearson') {       nbR[,gn]=residuals(fnbrF,'pearson') }
        cisNB[[gn]]$transcript=gn
        cisNB[[gn]]$lmarker=p
        cisNB[[gn]]$theta=theta.est
        cisNB[[gn]]$negbin.beta=as.numeric(coef(fnbrF)[msize])
        cisNB[[gn]]$negbin.se=as.numeric(fnbrF$se[msize])
        cisNB[[gn]]$LLik=as.numeric(logLik(fnbrF))
        cisNB[[gn]]$negbin.LRS=as.numeric(-2*(logLik(fnbrN)-logLik(fnbrF)))
        cisNB[[gn]]$negbin.p=pchisq( cisNB[[gn]]$negbin.LRS,1, lower.tail=F) 
         #-2*(nmodel-plr)
        cisNB[[gn]]$fmodelBs=coef(fnbrF)
        
        if(resids=='pearson') {       cisNB[[gn]]$pearsonResiduals=residuals(fnbrF,'pearson') }
        
        #print(cisNB[[gn]])
    }

   # nbR.var=colVars(nbR)
   # nbR=nbR[,-which(is.na(nbR.var))]
    return(cisNB) #, nbR=nbR))
}

# formula for calculating the negative binomial glm log likelihood 
nbLL=function (y, mu, theta)  {
    return( -sum ( (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta + y) ))
}
# fast calculation of negative binomial LL across the genome if theta is known 
domap=function(gn, ...) { 
       theta.est=thetas[gn]
       YY=Y[,gn]
       nbn=negative.binomial(theta=theta.est,link='log')
      
       fnbrN=fastglmPure(DM, YY, family=nbn)
       nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)

       LRS=rep(NA,ncol(Gsub))
       names(LRS)=colnames(Gsub)

       #initialize the last column of DM to be 1s
       XXn=cbind(DM,1)
       idx=1:ncol(Gsub)
       gcovn=(ncol(XXn))
       for(gx in idx) { 
           #colnames(Gsub)){
           XXn[,gcovn]=Gsub[,gx]
           fnbrF=fastglmPure(XXn, YY,  family=nbn)
           LRS[gx]=-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))
        }
       return(LRS)
}


cl <- makeCluster(36)
clusterEvalQ(cl, {
       library(fastglm)
       library(Matrix)
       library(MASS) 
       NULL  })


#cc=read.csv(paste0(reference.dir, 'cell_cycle.csv')) #read.csv(paste0(


LDprune=function(geno, m.granges, cor.cutoff=.999) {
    G=standardise2(vg)
    #prune markers 
    chr.ind=split(1:length(m.granges), as.character(seqnames(m.granges)))
    G.chr=lapply(chr.ind, function(x) G[,x])
    Gcor.chr=lapply(G.chr, function(x) crossprod(x)/(nrow(x)-1))
    fcm.chr=lapply(Gcor.chr, function(x) caret::findCorrelation(x,cutoff=cor.cutoff, verbose=F) )
    fcm=as.vector(unlist(mapply(function(x,y) x[y],  sapply(G.chr, colnames), fcm.chr)))
    m.to.keep=which(!colnames(vg)%in% fcm) #colnames(vg)[fcm])
    Gsub=vg[,m.to.keep]
    rm(Gcor.chr)
    m.to.keep=match(colnames(Gsub), colnames(vg))
    markerGRr=m.granges[m.to.keep]
    return(list(Gsub=Gsub, markerGRr=markerGRr))
}


#overwrite these functions
sgd.genes=getSGD_GeneIntervals(reference.dir)
getCisMarker=function(sgd.genes, m.granges, counts) {
    Yre=t(counts)
    markerGR=m.granges
    sgd.genes=sgd.genes[match(colnames(Yre), sgd.genes$gene_id),]
    tag=sgd.genes

    cisMarkers=data.frame(transcript=sgd.genes$gene_id, 
                          marker=m.granges$sname[nearest(tag,markerGR)],
                          stringsAsFactors=F)
    cisMarkers=na.omit(cisMarkers)
    return(cisMarkers)
    #indices of cis eQTL in 
    #rows are transcripts, columns are markers
    #cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
}

doCisNebula=function(gn,...){
        mmpN=cbind(mmp1, Gsub[,cSplit[[gn]]$marker[1]])
        yind=match(cSplit[[gn]]$transcript, rownames(counts))
        NBcis=nebula(counts[yind,], mmpN[,1], pred=mmpN)
        NBcis$summary$gene=rownames(counts)[yind]
        return(NBcis)
}


nUMI_thresh=10000
minInformativeCellsPerTranscript=128
#for(ee in c(2,4,7:12,15,16)) { #1:length(experiments)){
for(ee in c(4,7:10,11,12)) { #,15,16)) { #1:length(experiments)){
#for(ee in c(15,16)) { #1:length(experiments)){

    experiment=experiments[ee]
    print(ee)
    print(experiment)    
    #cross=crossL[[cdesig[ee]]]
    #parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')

    #read in transcript counts 
    counts=readRDS(paste0(results.dir, 'counts.RDS'))

    #read in geno informative counts
    g.counts=readRDS(paste0(results.dir, 'gcounts.RDS'))
    
    # read in genetic map and format for plotting
    gmapC=rbindlist(readRDS(paste0(results.dir, 'gmap.RDS')))
    gmap.subset=gmapC[match(rownames(g.counts[[1]]), paste0(gmapC$chrom, '_', gmapC$ppos)),]
    gmap.ss=split(gmap.subset, gmap.subset$chrom) 
   
    #get hmm genotype probs 
    vg=getSavedGenos(chroms, results.dir, type='genoprobs')
    m.granges=getMarkerGRanges(g.counts)
    m.granges$gcoord=gcoord.key[as.character(seqnames(m.granges))]+start(m.granges)
    m.granges$sname=colnames(vg)

    #add additional hard filter for umis per cell
    classification = readRDS(paste0(results.dir, 'cellFilter.RDS'))

    #add additional hard filter for umis per cell
    classification = classification & colSums(counts)<nUMI_thresh

    #matches data to cell.covariates
    counts=counts[,names(classification)[classification]] #cell.covariates$barcode]
    vg=vg[names(classification)[classification],] #cell.covariates$barcode,]
   
    pruned=LDprune(vg, m.granges)
    Gsub=pruned$Gsub
    markerGRr=pruned$markerGRr
    rm(pruned)
    
    # af diagnostic    
    png(file=paste0(results.dir, 'seg_diagnostic_af_plot.png'), width=1024,height=1024)
    par(mfrow=c(2,1))
    plot(m.granges$gcoord, colSums(vg)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='fraction of segs that are parent 2', 
         main=experiment, sub=paste(ncol(vg), 'markers'))
    abline(v=gcoord.key, lty=2, col='lightblue')
    abline(h=0.5, col='grey')
    #dev.off()
    #png(file=paste0(results.dir, 'seg_diagnostic_af_plot_pruned.png'), width=1024,height=512)
    plot(markerGRr$gcoord, colSums(Gsub)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='fraction of segs that are parent 2',
         main=experiment, sub=paste(ncol(Gsub), 'markers after pruning'))
    abline(v=gcoord.key, lty=2, col='lightblue')
    abline(h=0.5, col='grey')
    dev.off()
   
    rm(vg) 
    # relatedness between segs 
    tG=t(Gsub)
    tt=hclust(as.dist(1-(crossprod(scale(t(Gsub)))/(nrow(tG)-1))^2))
    png(file=paste0(results.dir, 'seg_diagnostic_corr_clust_plot.png'), width=1024,height=512)
    plot(tt, labels=F, main=experiment, ylab="1-r^2")
    dev.off()
    
    #count informative cells per transcript
    infoCells=rowSums(counts>0)
    expressed.transcripts=names(infoCells)[(infoCells>minInformativeCellsPerTranscript)] #names(sum(infoCells>128)) #[infoCells>128,]

    #select only the transcripts expressed in enough cells 
    Yr=t(counts[expressed.transcripts,])
    
    cisMarkers=getCisMarker(sgd.genes[sgd.genes$gene_id %in% rownames(counts),], markerGRr, t(Yr)) #counts)
    #cisMarkers=cisMarkers[cisMarkers$transcript %in% expressed.transcripts,]
    #if no subsetting this is not necessary
    #cMsubset=cisMarkers[which(cisMarkers$transcript %in% colnames(Yr)),]

    transcript.features=data.frame(gene=colnames(Yr), non.zero.cells=colSums(Yr>1), tot.expression=colSums(Yr))
    barcode.features=data.frame(barcode=colnames(counts), nUMI=colSums(counts))

    mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts))))
    cSplit=split(cisMarkers, cisMarkers$marker)
    #dispatch is slower than execute 
    clusterExport(cl, varlist=c("mmp1", "Gsub", "cSplit", "counts","doCisNebula", "nebula"))
    cisNB=parLapply(cl, names(cSplit), doCisNebula)
    names(cisNB)=names(cSplit)

    #cisNB=mclapply(cSplit, function(x,...) {
    #    mmpN=cbind(mmp1, Gsub[,x$marker[1]])
    #    yind=match(x$transcript, rownames(counts))
    #    NBcis=nebula(counts[yind,], mmpN[,1], pred=mmpN)
    #    NBcis$summary$gene=rownames(counts)[yind]
    #    return(NBcis)
    #    },mmp1,cisMarkers,Gsub,counts, mc.cores=24)
    saveRDS(cisNB, file=paste0(results.dir, 'cisNB2.RDS'))
    cis.ps=as.vector(unlist(sapply(cisNB, function(x) x$summary$p_)))
    names(cis.ps)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))

    png(file=paste0(results.dir, 'cis_p_hist.png'), width=512,height=512)
    hist(cis.ps,breaks=50, main=experiment, sub=sum(p.adjust(cis.ps,'fdr')<.05))
    dev.off()

    #here we define covariates 
    thetas=as.vector(1/unlist(sapply(cisNB, function(x) x$overdispersion$Cell)))
    names(thetas)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))

    nullNB=nebula(counts[expressed.transcripts,], mmp1[,1], pred=mmp1)
    dispersion.df=data.frame(gene=nullNB$summary$gene, theta.null=1/nullNB$overdispersion[,2])
    cisModel.df=data.frame(gene=names(thetas),theta.cis=thetas)
    dispersion.df=left_join(dispersion.df, cisModel.df, by='gene')
    saveRDS(dispersion.df, file=paste0(results.dir, 'dispersions.RDS'))

    #could switch between thetas here 
    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene

    print('Calculating LODs')
    #Yks=Matrix(Yr, sparse=T)
    clusterExport(cl, varlist=c("thetas", "Yr", "Gsub", "mmp1", "nbLL", "domap"))
    clusterEvalQ(cl, { Y=Yr;    DM=mmp1;   return(NULL);})
    LOD=do.call('rbind', parLapply(cl, names(thetas), domap) )
    LOD[is.na(LOD)]=0
    LOD=LOD/(2*log(10))
    rownames(LOD)=names(thetas)
    saveRDS(LOD, file=paste0(results.dir, 'LOD_NB3.RDS'))
    msig=which(LOD>5, arr.ind=T)
    png(file=paste0(results.dir, 'seg_diagnostic_plot_3.png'), width=768,height=768)
    plot(msig[,2], msig[,1])
    dev.off()
}



    #insert code speedup for dispersion estimate here
    L=irlba(Gsub,10)
    colnames(L$u)=paste0('u',1:ncol(L$u))
    mmpN=model.matrix(lm(Yr[,1]~log(colSums(counts))))
    mmpN=cbind(mmpN, L$u)

    pcNB=nebula(counts[expressed.transcripts,], mmpN[,1], pred=mmpN)
    pcModel.df=data.frame(gene=pcNB$summary$gene,theta=1/pcNB$overdispersion[,2])
    nullModel.df=data.frame(gene=pcNBnull$summary$gene,theta=1/pcNBnull$overdispersion[,2])

    cisModel.df=data.frame(gene=names(thetas),theta=thetas)

    test=left_join(cisModel.df, pcModel.df, by='gene', suffix=c('.cis', '.pc'))
    test=left_join(test, nullModel.df, by='gene', suffix=c('','.null'))
    theta_v_expression=left_join(transcript.features, test, by='gene')
    
    library(viridis)
    library(ggpubr)
    a=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta.cis),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()
    b=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta.pc),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()

  
   cc=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()

    ggarrange(a,b,cc)
plot(log(test$theta.x), log(test$theta.y))

    saveRDS(pcNB, file=paste0(results.dir, 'pcNB3.RDS'))

    pcNBnull=nebula(counts[expressed.transcripts,], mmp1[,1], pred=mmp1)


    test=cv.glmnet(mmpN, Yr[,4], family='poisson', alpha=1)

    
    print('Calculating overdispersions') 
    pcNB=nebula(counts[expressed.transcripts,], mmpN[,1], pred=mmpN)
    saveRDS(pcNB, file=paste0(results.dir, 'pcNB3.RDS'))
    #pcNB=readRDS(paste0(results.dir, 'pcNB2.RDS'))
    thetas=1/pcNB$overdispersion$Cell  
    names(thetas)=pcNB$summary$gene



    theta_v_expression=left_join(pcModel.df, transcript.features, by='gene')
    ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()
  









       # LODo=readRDS(paste0(results.dir, 'LOD_NB2.RDS'))


}









match(rownames(LOD), sgd.genes$gene_id)


    resids=sapply(cisNB, function(x) x$pearsonResiduals)
    rownames(resids)=rownames(Yr)

    cc.resids=resids[,colnames(resids) %in% cc$ORF]
    um=uwot::umap(cc.resids, n_neighbors=30, metric="cosine", min_dist=0.3, n_threads=36)
    df=data.frame(u1=um[,1], u2=um[,2],data.matrix(Yr))
    ggplot(df,aes(x=u1,y=u2, col=log10(YBL003C)+1))+scale_color_viridis()+geom_point()

    rawC=standardise2(data.matrix(log(Yr[,colnames(Yr) %in% cc$ORF]+1)/log(colSums(counts))))
    um2=uwot::umap(rawC, n_neighbors=30, metric="cosine", min_dist=0.3, n_threads=36)




    df=data.frame(u1=um2[,1], u2=um2[,2],data.matrix(Yr))
    ggplot(df,aes(x=u1,y=u2, col=log10(YBL003C)+1))+scale_color_viridis()+geom_point()





 #   cell.cycle.classification.file=paste0(base.dir, cell.cycle.classification.dir, cell.cycle.classification.names[ee], '.cluster.assignments.txt')
 #   cell.cycle.annot.file=paste0(base.dir, cell.cycle.classification.dir, cell.cycle.classification.names[ee], '.cluster.annotations.txt')
    # get cell cycle covariates
#    cell.covariates=getCCinfo(cell.cycle.classification.file,cell.cycle.annot.file)
    #note n gets smaller here due to cells missing classification
 #   cell.covariates=cell.covariates[!is.na(cell.covariates$cid),]
 #   print(experiment)
 #   print(nrow(cell.covariates))
 #   saveRDS(cell.covariates, paste0(results.dir, 'cell_covariates.RDS'))
  
    
    # par(xaxs='i', yaxs='i')
  #  plot(m.granges$gcoord, colSums(pg[rownames(pg)%in%rownames(vg),])/2727, ylim=c(0,1))
  #  plot(m.granges$gcoord, colSums(vg[classification,]-1)/sum(classification), ylim=c(-1,1), xlab='mpos', ylab='af')
  #  abline(v=gcoord.key, lty=2, col='lightblue')
  #  cis.marker=getCisMarker(sgd.genes, m.granges, counts)
 #vgV=getSavedGenos(chroms, results.dir, type='viterbi')
    #par(xaxs='i', yaxs='i')
    #plot(m.granges$gcoord, colSums(pg[rownames(pg)%in%rownames(vg),])/2727, ylim=c(0,1))
    #plot(m.granges$gcoord, colSums(vgV[classification,]-1)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='af')
    #abline(v=gcoord.key, lty=2, col='lightblue')

     # assuming classification matches experiment, subset and propagate order across phenos and genos 

 #cisNB=fitNBCisModel(cMsubset, mmp1, Yr, Gsub, resids='pearson')          #, resids='deviance')
    #saveRDS(cisNB, file=paste0(results.dir, 'cisNB.RDS'))
    #thetas=sapply(cisNB, function(x) as.numeric(x$theta))
    #intercepts=sapply(cisNB, function(x) as.numeric( x$fmodelBs[1]) )
    #plot(1/comp_models$theta.x, 1/comp_models$theta.y, main='1/theta', xlim=c(0,25), ylim=c(0,25), xlab='cis model', ylab='10 pc model')
    
    #cisModel=readRDS(paste0(results.dir, 'cisNB.RDS'))
    #pcModel=readRDS(paste0(results.dir, 'pcNB2.RDS'))
    #thetas=sapply(cisModel, function(x) as.numeric(x$theta))
    #cisModel.df=data.frame(gene=names(thetas), theta=thetas)
    #pcModel.df=data.frame(gene=pcModel$summary$gene,theta=1/pcModel$overdispersion[,2])
    #comp_models=left_join(cisModel.df, pcModel.df, by='gene')











    ## optional get ASE info
    # ase.counts=getASEcounts(g.counts, sgd.genes, m.granges, 16)
    # tot.geno.inf=ase.counts[[1]]+ase.counts[[2]] 

    # taf=colSums(ase.counts$ref)/(colSums(ase.counts$ref+ase.counts$alt))
    # caf=log2(colSums(ase.counts$ref+ase.counts$alt))
    # plot(sgd.genes$gcoord[match(names(taf), sgd.genes$ID)],taf, col=(caf<5)+1)
    # abline(v=gcoord.key, lty=2, col='lightblue')
    # hist(log2(colMeans(tot.geno.inf)), xlab='log2(avg genotype informative counts per transcript per cell)')


    #Gcor=crossprod(G)/(nrow(G)-1)
    #rm(G)
    # visualize auto-correlation. e.g. correlation decay from given marker to a marker with an index 100 less
    #v=Gcor[row(Gcor)==(col(Gcor)-100)] ... eek
    # prune correlated markers  (approx every 5cm)
    # previous code transformed genotypes as log(g/(1-g)) for a beta-regression-like setup
    # ~30 or markers were lost by div by 0 errors, fixed 2/14
    #fcm=caret::findCorrelation(Gcor,cutoff=.999, verbose=F)














 
    ## get genotyypes
    #vg=getSavedGenos(chroms, results.dir, type='viterbi')
    #saveRDS(vg, paste0(results.dir, 'vit_gmatrix.RDS'))

    #pg=getSavedGenos(chroms, results.dir, type='genoprobs')
    
        # antiquated, read in r/qtl binned formatted genotypes   
    #rQTL.coded.file=paste0(results.dir, 'rQTLcoded.RDS')
    #rQTL.coded=readRDS(rQTL.coded.file)


















    #plot(m.granges$gcoord, colSums(vg-1)/nrow(vg), ylim=c(0,1))
    #abline(v=gcoord.key)

    cell.covariates$total.counts=colSums(counts)
    
    Yr=t(as.matrix(counts))
    Yr.var=colVars(Yr, parallel=T)
    Yr=Yr[,Yr.var>0]
    tcounts=Yr>0
    tcounts=colSums(tcounts)
    Yr=Yr[,tcounts>20]

    countsR=counts[colnames(Yr),]

    #.05 230
    mm=model.matrix(counts[1,]~log2(total.counts)+cid,data=cell.covariates) 
    #Yre2=matrix(NA, nrow(Yr), ncol(Yr))
    #colnames(Yre2)=colnames(Yr)
    #rownames(Yre2)=rownames(Yr)
    #for(g in 1:ncol(Yr)){
    #    print(g)
    #    Yre2[,g]=scale(residuals(glm(Yr[,g]~offset(log(total.counts))+cid,data=cell.covariates, family=poisson()),'pearson' ))
    #}
    
    #mm0=model.matrix(counts[1,]~log2(total.counts),data=cell.covariates) 
    #BOLS=lm.fit(mm0, log2(Yr+1))
    #Yr0=scale(residuals(BOLS))
    #rPC=hd.eigen(t(Yr0),vectors=T, center=F, scale=F, k=2)
    #plot(rPC$vectors[,1], rPC$vectors[,2])

    BOLS=lm.fit(mm, log2(Yr+1))
    Yresid=residuals(BOLS)
    rownames(Yresid)=rownames(Yr)
    colnames(Yresid)=colnames(Yr)
    rm(BOLS)

    Yre=standardise2(Yresid)
    rownames(Yre)=rownames(Yresid)
    colnames(Yre)=colnames(Yresid)
    #------------------------------------------------------------------------------------------------------------

    #Yre=Yre2
    assignme=model.matrix(counts[1,]~cell.covariates$cid-1)
    cgen=scale(residuals(lm(assignme~log2(cell.covariates$total.counts))))

    rcc=crossprod(cgen,G)/(nrow(Yre)-1)
    LODcc=-nrow(Yre)*log(1-rcc^2)/(2*log(10))
    #saveRDS(LODcc, file=paste0(results.dir, 'lodCC.RDS'))
    par(mfrow=c(7,1), xaxs='i')
    for(i in 1:7) {
    plot(m.granges$gcoord, LODcc[i,], main=rownames(rcc)[i], ylab='LOD')
    abline(v=gcoord.key, lty=2, col='lightblue')
    }


   # 1D scan speedup 
    r=crossprod(Yre,G)/(nrow(Yre)-1)
    LOD=-nrow(Yre)*log(1-r^2)/(2*log(10))

    #pool for permutations
    wgFDR= getFDRfx(r, Yre, G, nperm=5)
    peaks1D=get1Dpeaks(chroms,r,LOD,wgFDR, m.granges, sgd.genes)
    saveRDS(peaks1D, file=paste0(results.dir, 'peaks1D.RDS'))
    plot2D(peaks1D, .1, gcoord.key=gcoord.key, experiment=experiment)
   
   
    #multivariate scan for hotspots 
    mLOD=list()
    for(cc in chroms) {
        print(cc)
        ng=sgd.genes[seqnames(sgd.genes)!=cc,]
        YreS=Yre[,(colnames(Yre) %in% ng$Name)]
        moi=m.granges$sname[as.character(seqnames(m.granges))==cc]
        #test=hd.eigen(t(Yre), center=F, scale=F,vectors=T)
        pc.to.retain=40
        yreduced=hd.eigen(t(YreS), center=F, scale=F, k=pc.to.retain, vectors=T)
        testF=mvn.scanone(G[,moi], yreduced$vectors)
        testN=determinant(crossprod(yreduced$vectors), logarithm=T)$modulus
        mLODv=(nrow(yreduced$vectors)/2)*(testN-testF)/(2*log(10))
        mLOD[[cc]]=mLODv
    }
    mLOD=do.call('c', mLOD)    
    peaksMV=getMVpeaks(chroms, r, LOD, mLOD, wgFDR, m.granges, sgd.genes)
    #saveRDS(peaksMV, file=paste0(results.dir, 'peaksMV.RDS'))


    ccPeaks=list()
    # split by cell -cycle
    cell.covs=split(cell.covariates, cell.covariates$cid)
    #Yr and vg are counts and genos respectively
    for(covn in names(cell.covs)){
            print(covn)
          ccdf=cell.covs[[covn]]
          Yrs=Yr[ccdf$barcode,]
          tcounts=Yrs>0
          tcounts=colSums(tcounts)
          Yrs=Yrs[,tcounts>20]
          
          mms=model.matrix(Yrs[,1]~log2(total.counts),data=ccdf) 
          BOLSs=lm.fit(mms, log2(Yrs+1))
          Yresids=residuals(BOLSs)
          rm(BOLSs)
          rownames(Yresids)=rownames(Yrs)
          colnames(Yresids)=colnames(Yrs)
          Yres=standardise2(Yresids)
          rownames(Yres)=rownames(Yresids)
          colnames(Yres)=colnames(Yresids)
        
          vgs=vg[ccdf$barcode,]
          Gs=standardise2(vgs)
          rs=crossprod(Yres,Gs)/(nrow(Yres)-1)
          LODs=-nrow(Yres)*log(1-rs^2)/(2*log(10))

         #pool for permutations
         wgFDRs= getFDRfx(rs, Yres, Gs, nperm=5)
         peaks1Ds=get1Dpeaks(chroms,rs,LODs,wgFDRs, m.granges, sgd.genes)
         ccPeaks[[covn]]=peaks1Ds
    }
    saveRDS(ccPeaks, file=paste0(results.dir, 'CCpeaks.RDS'))



}
texp=colSums(Yr)
llr=apply(ll,1,range)
llm=apply(ll,1,mean)
plot(log10(texp), max.obsLOD,col='red', ylim=c(0, 0.2),xlim=c(1,6.5),
     xlab='log10(total expression per transcript)',
     ylab='abs(r)'
    )
points(log10(texp), llm, ylim=c(0.02,0.06))
segments(log10(texp),llr[1,],log10(texp),llr[2,], col="#00000066")
points(log10(texp), max.obsLOD,col='red', ylim=c(0, 0.2))
#points(log10(texp), llm,  col="#00000022")

qm=rbindlist(ccPeaks, idcol='cell_cycle')
qm=qm[qm$FDR<.1,]
ggplot(qm,aes(x=mgcoord,y=tgcoord, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x =mgcoordL, y = tgcoord, xend =mgcoordR, yend = tgcoord)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    facet_wrap(~cell_cycle)+
    geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/100919/yeast_CC_QTL.png')


countsR=as.matrix(countsR)

p1dr=peaks1D[peaks1D$FDR<.1,]
p1dr$transcript=droplevels(p1dr$transcript)
gqtl=split(p1dr, p1dr$transcript)

cell.covariates$cid=relevel(cell.covariates$cid, "G1")

cc.fits=foreach(gene=names(gqtl)) %do% {
    # cell cycle 
    tryCatch({
    #gene=names(gqtl)[3]

    y=countsR[gene,]
    XX=G[,as.character(gqtl[[gene]]$peak.marker)]
    #XX=matrix(XX)
    #ZXX=Z%*%XX
    #colnames(XX)=as.character(gqtl[[gene]]$peak.marker) #sb[[gene]]$fscan.marker
    yp=log2(y+1)

    obs=data.frame(y=y, yp=yp, cell_total=cell.covariates$total.counts, Cid=cell.covariates$cid)

    #m2=model.matrix(y~cell.covariates$cid-1)
    #mm0=lm(yp~log2(cell_total)+XX,data=obs)
    #mm0=lm(yp~log2(cell_total)+Cid+XX,data=obs)

    #given this is multithreaded, do not use dopar!
    mm=lm(yp~log2(cell_total)+XX*Cid,data=obs)
    #mm=glm.nb(y~offset(log(cell_total))+XX*XX*Cid,data=obs)
    #mm=glm.nb(y~(log(cell_total))+XX*Cid,data=obs)

    clist=tidy(mm)
    if(is.matrix(XX)==FALSE) {
      clist$term=gsub('XX',   paste0('XX', as.character(gqtl[[gene]]$peak.marker)), clist$term)
    }
    print(gene)
    print(data.frame(clist))
     return(clist)
    },error=function(e) {
        return(NULL)
    })
}
names(cc.fits)=names(gqtl)
saveRDS(cc.fits, file = paste0(results.dir, 'CCfits.RDS'))
#save(cc.fits, file='/data/single_cell_eQTL/yeast/cc.fits.RData')
# for regular lm
cc.ps=sapply(cc.fits, function(x) {
                 y=x$p.value[grep(':', x$term)] 
                 names(y)=x$term[grep(':', x$term)]
                 return(y)
    } )
hist(unlist(cc.ps),breaks=350, xlab='p', main='(QTL x Cell_Cycle) p-values')

rcc.ps=lapply(cc.ps,function(n) {
           ss=strsplit(names(n),':')
           qm=sapply(ss, function(x) gsub('XX', '', x[1]))
           cm=sapply(ss, function(x) x[2])
           data.frame(peak.marker=qm,cell_cycle=cm,p=as.vector(n),stringsAsFactors=F)
    })
rcc.ps=rbindlist(rcc.ps, idcol='transcript')
rcc.ps$FDR=qvalue(rcc.ps$p)$qvalue

    rcc.ps$mchr=as.character(seqnames(m.granges))[match(rcc.ps$peak.marker, m.granges$sname)]
    rcc.ps$mpos=(start(m.granges))[match(rcc.ps$peak.marker, m.granges$sname)]
    rcc.ps$mgcoord =m.granges$gcoord[match(rcc.ps$peak.marker, m.granges$sname)]
    #rcc.ps$mgcoordL=m.granges$gcoord[match(rcc.ps$CI.l, m.granges$sname)]
    #rcc.ps$mgcoordR=m.granges$gcoord[match(rcc.ps$CI.r, m.granges$sname)]
    #rcc.ps$mposL   =(start(m.granges))[match(rcc.ps$CI.l, m.granges$sname)]
    #rcc.ps$mposR   =(start(m.granges))[match(rcc.ps$CI.r, m.granges$sname)]

    rcc.ps$tchr=as.character(seqnames(sgd.genes))[match(rcc.ps$transcript, sgd.genes$Name)]
    rcc.ps$tpos=start(sgd.genes)[match(rcc.ps$transcript, sgd.genes$Name)]
    rcc.ps$tgcoord=sgd.genes$gcoord[match(rcc.ps$transcript, sgd.genes$Name)]
saveRDS(rcc.ps, file=paste0(results.dir, 'CCpeaksxQTL_int.RDS'))

qm2=rcc.ps[rcc.ps$FDR<.1,]

ggplot(qm2,aes(x=mgcoord,y=tgcoord, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    #geom_segment(aes(x =mgcoordL, y = tgcoord, xend =mgcoordR, yend = tgcoord)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    facet_wrap(~cell_cycle)+
    geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())           

sq=split(qm2, qm2$transcript)
sq=sq[which(sapply(sq,nrow)>1)]
sqf=lapply(sq,function(x) x[x$mchr!='chrIV',])
t1=cc.fits[names(which(sapply(sqf,nrow)>1))]

#svF=irlba(G, nv=150)
#svFGG=irlba(G*G, nv=150)

disp.fits=foreach(gene=names(gqtl)) %dopar% {

 tryCatch({
    #gene=rownames(counts)[r]
    #  dispersion test 
    #gene='YHB1'
    #gene='YJR009C'
    print(gene)
    #y=raw_data[r,]
    y=countsR[gene,]
    XX=G[,as.character(gqtl[[gene]]$peak.marker)]
    yp=log2(y+1)
    obs=data.frame(y=y, yp=yp, cell_total=cell.covariates$total.counts, Cid=cell.covariates$cid)
    #test=model.matrix(y~(XX[,1]+XX[,2]+XX[,3]+XX[,4])^2-1)
    bm=glmmTMB(y~offset(log(cell_total))+Cid+XX, family=nbinom2,data=obs)
    #bm2=glmmTMB(y~offset(log(cell_total))+Cid+test, family=nbinom2,data=obs)

    # bm=glmmTMB(y~(log(cell_total))+Cid+XX, family=nbinom2,data=obs)
    #bm2=update(bm, dispformula=~Cid+XX)


   if(is.null(dim(XX)[2])) {XX=as.matrix(XX)}
   dlist=list()
   for(i in 1:ncol(XX)){
            bm1=update(bm, dispformula=~XX[,i])
            ps=as.numeric(pchisq(2*(logLik(bm1)-logLik(bm)),1,lower.tail=F))
            #nvec=c(summary(bm1)$coefficients$cond[i+4,], sigma(bm1), ps)
            nvec=c(summary(bm1)$coefficients$cond[i+4,], sigma(bm1),ps,summary(bm1)$coefficients$disp[2,])
            names(nvec)[5]='theta_int'
            names(nvec)[6]='theta_coef'
            names(nvec)[7]='lrt_p'
            names(nvec)[8:11]=paste0('theta_coef_',names(nvec)[8:11]) 
            dlist[[as.character(gqtl[[gene]]$peak.marker)[i]]]=nvec
        }

    print(dlist)
    return(dlist)

    },error=function(e) {
        dlist=list()
        nvec=rep(NA,11)
        names(nvec)=c("Estimate",  "Std. Error",            "z value",               "Pr(>|z|)"  ,           
         "theta_int"  ,           "theta_coef"      ,      "lrt_p"   ,              "theta_coef_Estimate" , 
         "theta_coef_Std. Error", "theta_coef_z value"  ,  "theta_coef_Pr(>|z|)")  
        dlist[['fail']]=nvec
        return(dlist)
    })
}
names(disp.fits)=names(gqtl)
hist(unlist(sapply(disp.fits, function(x) sapply(x, function(y)y[4]))),breaks=100)
hist(unlist(sapply(disp.fits, function(x) sapply(x, function(y)y[7]))),breaks=100)


dp3df=rbindlist(lapply(disp.fits, function(x) data.frame(t(x[[1]]))),idcol='gene')

df.table=rbindlist(disp.fits,idcol='transcript',fill=T)

df.table=rbindlist(disp.fits,idcol='transcript')




    bm=glmmTMB(y~offset(log(cell_total))+Cid*XX, family=nbinom2,data=obs)
    clist=tidy(bm)
    bm2=glm.nb(y~offset(log(cell_total))+Cid*XX,data=obs)
    
    #counts[1,]~log2(total.counts)+cid,data=cell.covariates) 
    



h=t2$cPeaks[t2$cPeaks$chr=='chrVIII',]
par(xaxs='i')
plot(m.granges$gcoord, mLOD, type='l')
abline(v=gcoord.key, lty=2, col='lightblue')

     
# chrVIII is gpa1
# ylim=c(0,30), type='l')  # /4.6))


plot2D=function(peaks1D, FDR.thresh=.05, gcoord.key=gcoord.key, experiment=experiment) {
    plot.thresh=FDR.thresh
    par(xaxs='i', yaxs='i')
    plot(peaks1D$mgcoord[peaks1D$FDR<plot.thresh], peaks1D$tgcoord[peaks1D$FDR<plot.thresh],
         ylab='transcript position',
         xlab='QTL position', main=experiment, xaxt='n', yaxt='n', pch=20)
    axis(1, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    axis(2, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    abline(v=gcoord.key, lty=2, col='grey')
    abline(h=gcoord.key, lty=2, col='grey')
}







peaks1D=readRDS('/data/single_cell_eQTL/yeast/results/2_2444_44-/peaks1D.RDS')


 


















    # a chromosome of interest for plotting 
    coii='chrVII'

    pp.file=paste0(results.dir, coii,'.RDS')
    post.prob=readRDS(pp.file)
    
    v.file=paste0(results.dir, coii, '_viterbi.RDS')
    vit= readRDS(v.file)

    for(i in 100:200) {
        diagnoseHMM(i,chrom=coii, gmap.ss, g.counts[[1]], g.counts[[2]], post.prob, rQTL.coded, viterbiPath=vit, classification=classification)
        readline()
    }

}
















    # recurring.ref=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==1,na.rm=T))
    # recurring.alt=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==2,na.rm=T))
    # raf=recurring.ref/(recurring.alt+recurring.ref)
    # taf=recurring.ref+recurring.alt
    # approx rate of hets vs homzyg = 0.0075
    # 29234/(2016111 +1881594)



    #saveRDS(classifications, file = paste0(base.dir, 'is_haploid_classifcation.RDS'))
#################################################

 
# usually we'll want to parallelize this further

#v=viterbi(cross2, error_prob=.01, ,lowmem=F,quiet=F, cores=16)
#dv=do.call('cbind', v)
#by.freq=apply(dv,2, function(x) sum(x==1)/length(x))
#dv.prune=dv[,by.freq<.7 & by.freq>.3]
#G=dv.prune

G2=do.call('cbind', sapply(gps, function(x) x[,2,]))

YJMpS=rowSums(G2)


Ys=t(as.matrix(counts))
Ys=Ys[-het.cells,]

total.counts=rowSums(Ys)

mm=model.matrix(Ys[,1]~log2(total.counts))
BOLS=lm.fit(mm, log2(Ys+1), intercept=F)
Yresid=residuals(BOLS)
rownames(Yresid)=rownames(Ys)
colnames(Yresid)=colnames(Ys)

Yre=standardise(Yresid)
#Yre=standardise(Rfast::colRanks(Yresid))
rownames(Yre)=rownames(Yresid)
colnames(Yre)=colnames(Yresid)
nac=apply(Yre, 2, function(x) sum(is.na(x)))
Yre=Yre[,nac==0]


Gr=standardise(G2)
colnames(Gr)=colnames(G2)
rownames(Gr)=rownames(Yre)

r=crossprod(Yre,Gr)/(nrow(Yre)-1)
LOD=(-nrow(Yre)*log(1-r^2))/(2*log(10))
sigL=which(LOD>8, arr.ind=T)
plot(sigL[,2], sigL[,1], pch=21)

getFDRfx=function(r, Yre,Gr,nperm=5, vint=seq(0.001,.15,.001)){
    
    max.obsLOD=apply(abs(r),1,max)

    permR=replicate(nperm, {
                    nind=sample(1:nrow(Yre))
                    crossprod(Yre[nind,],Gr)/(nrow(Yre)-1)
    })

    ll=apply(abs(permR),3, function(x) apply(x,1,max))
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = vint

    # expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(apply(ll,2,function(x) sum(x>thresh))))
                 })
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    

    pFDR = rev(cummax(rev(pFDR)))
    fdrFX=approxfun(pFDR, vint, ties=mean)
    rtoFDRfx=approxfun(vint,pFDR, ties=mean)

    return(list(fdrFX=fdrFX,rtoFDRfx=rtoFDRfx))
}
#getPeaks=function(r, thresh){
#    cind=gsub('_.*','', colnames(r))
#    cmaxind=matrix(0, nrow(r), 16)
#    rownames(cmaxind)=rownames(r)
#    cmax=cmaxind
#    for(i in 1:nrow(r)) {
#        #print(i)
#        x=r[i,]
#      y=split(x, cind);
#     cmaxind[i,]=as.vector(sapply(y, function(h) names(which.max(abs(h)))))
#     cmax[i,]=x[cmaxind[i,]]
#    }
#    tmax=rep(rownames(cmaxind),16)
#    sig.hits=which(abs(cmax)>thresh) #,arr.ind=T)
#
#    qtl.p=cmaxind[sig.hits]
#    qtl.t=tmax[sig.hits]
#    return(data.frame(transcript=qtl.t,marker=qtl.p, stat=cmax[sig.hits]))
#}

#Lthresh=(nrow(Yre)*log(1-ff(.05))^2)/(2*log(10))
chromFDR=sapply(chroms, function(x) {
                 print(x)
                 Gs=Gr[,grep(paste0('^',x,'_'), colnames(Gr))]
                 r=crossprod(Yre,Gs)/(nrow(Yre)-1)
                 getFDRfx(r, Yre, Gs, nperm=5)
                 })


cPeaks=list()
for(cc in chroms) {
    print(cc)
    rS=r[,grep(paste0('^',cc,'_'), colnames(r))]
    mstat=apply(rS,1,max)
    mstat.pos=apply(rS,1,which.max)
    cPeaks[[cc]]=data.frame(
    transcript=rownames(rS),
    marker=colnames(rS)[mstat.pos],
    stat=mstat,
    FDR=chromFDR[,cc][[2]](mstat)
    )
    cPeaks[[cc]]$FDR[is.na(cPeaks[[cc]]$FDR)]=1
}


cP=do.call('rbind', cPeaks)
plot(match(cP$marker[cP$FDR<.1], colnames(Gr)),
     match(cP$transcript[cP$FDR<.1], colnames(Yre) ), 
             xlab='marker index', ylab='transcript index', main='joint analysis FDR < 10% corrected counts')

#viterbi = 894 at fdr<10%
