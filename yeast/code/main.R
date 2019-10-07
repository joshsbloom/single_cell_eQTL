library(monocle3)
library(Matrix)
library(rtracklayer)
library(data.table)
library(qtl)
library(qtl2)
library(parallel)
library(Rfast)
#library(mclust)
library(seqinr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(foreach)
library(doMC)
# for parallelizing the HMM
ncores=64
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
sgd=as.data.frame(import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff')))
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

# for now load all crosses--------------------------------------------------
load(paste0(reference.dir, 'cross.list.RData'))
load(paste0(reference.dir, 'parents.list.RData'))
#reorganize
crossL=list('A'=cross.list[['A']], 'B'=cross.list[['B']])
parentsL=list('A'=parents.list[['A']], 'B'=parents.list[['B']])
rm(cross.list)
rm(parents.list)

# need genetic map
genetic.maps=lapply(crossL, extractGeneticMapDf)
#----------------------------------------------------------------------------

# set some variables---------------------------------------------------------
base.dir='/data/single_cell_eQTL/yeast/'
cranger.dir='filtered_feature_bc_matrix/'
experiments=c('1_2444_44_1-2', '2_2444_44-', '3_ByxRM_51-_1-2', 
              '4_ByxRM_51-', '0_BYxRM_480MatA_1', '0_BYxRM_480MatA_2')
cdesig=c('B','B', 'A', 'A', 'A', 'A')
#by visual inspection (might want to revist this classification)
het.thresholds=c(2^6.1, 2^7.5, 2^7, 2^7.4, 2^6.5,2^6.5)
het.across.cells.threshold=.10

data.dirs=paste0(base.dir, experiments, '/')
classifications=list()
#-----------------------------------------------------------------------------

# iterate over each experiment -----------------------------------------------
for(ee in 1:length(experiments)){
    experiment=experiments[ee]
    cross=crossL[[cdesig[ee]]]
    parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    gmap=genetic.maps[[cdesig[ee]]]
    gmap.s=gmap; colnames(gmap.s)[2]='map'
    gmap.s=split(gmap.s, gmap.s$chrom)
    gmap.s=gmap.s[chroms]
    gmap.s=jitterGmap(gmap.s)

    results.dir=paste0(base.dir, 'results/', experiment, '/')
    dir.create(results.dir)

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
    #-------------------------------------------------------------

    # counts per parental allele (ref.counts and alt.counts, names not accurate) (fixes phasing) --------
    # also filters out counts from heavily distorted variants 
    g.counts=getGenoInformativeCounts(data.dir, cranger.dir, gmap,parents,3)
    #-----------------------------------------------------------------------------------------------------
  
    # can filter here for only sites that have informative variants in this data set,
    # but this may make comparison across datasets difficult in the future, hold on this for now
    rQTL.coded=encodeForRQTL(g.counts)
    
    het.count=apply(rQTL.coded, 2, function(x) sum(x==0,na.rm=T))
    tot.count=colSums(g.counts[[1]]+g.counts[[2]]) #ref.counts+alt.counts)
    #plot(log2(tot.count),log2(het.count))  # manual inspection of total counts vs sites classified as hets
    classifications[[experiment]]=het.count<het.thresholds[ee]
    saveRDS( classifications[[experiment]], file=paste0(results.dir, 'cellFilter.RDS'))

    ##png(file=paste0(data.dir, 'plot_classify_haploids.png'), width=1024, height=1024)
    ## add information about number of cells classified as haploid vs not
    #plot(log2(tot.count),log2(het.count), col=  classifications[[experiment]]+1,
    #     xlab='log2(total geno informative reads)', ylab='log2(total heterozygous sites)', 
    #     main=experiment, sub=paste('total = ', length(het.count),  ' haploid =', sum(classifications[[experiment]])))
    ##dev.off()

    het.cells=!classifications[[experiment]]
    recurring.hets=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==0,na.rm=T)) > (het.across.cells.threshold*sum(!het.cells))
    
    #per suggestion of James, run HMM on everything to start
    cross2=buildCross2(rQTL.coded, gmap, 'riself', het.sites.remove=T, het.cells.remove=F,
                       het.cells=het.cells, recurring.hets=recurring.hets)

    rqtl.genoprobs=calc_genoprob(cross2, error_prob=.005, lowmem=F, quiet=F, cores=16)
    rqtl.genoprobs=do.call('cbind', sapply(rqtl.genoprobs, function(x) x[,2,]))
    saveRDS(rqtl.genoprobs, file=paste0(results.dir, 'rqtl_genoprobs.RDS'))
    
    emissionProbs=estimateEmissionProbs(g.counts[[1]],  g.counts[[2]], error.rate=.002, 
                                recurring.het.sites.remove=T, 
                                recurring.hets=recurring.hets)
    
    runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='genoprob') #, n.indiv=1000)
    runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='viterbi') #, n.indiv=1000)

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
