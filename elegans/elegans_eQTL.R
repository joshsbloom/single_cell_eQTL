library(qtl)
library(qtl2)
library(qtl2convert)
library(data.table)
library(Rfast)

source('/data/single_cell_eQTL/elegans/elegans_eQTL_fx.R')
#aggregate phenotype data ------------------------------------------------------------------------
data.dir='/data/single_cell_eQTL/elegans/PE6Miseq/'
geneAnnot=readRDS(paste0(data.dir,'GeneAnnotationsPE6MiseqLibrary.rds'))
clusterID=readRDS(paste0(data.dir,'ClusterIdentityPE6MiseqLibrary.rds'))
counts=readRDS(paste0(data.dir, 'CountMatrixPE6MiseqLibrary.rds'))

# dcounts is a n X m matrix of umi counts per cell, rownames are cell_ids, 
# colnames are transcripts (ideally sorted by genomic position)
dcounts=data.matrix(t(as.matrix(counts)))

#total number of umi counts per cell
cell.counts=rowSums(dcounts)
#--------------------------------------------------------------------------------------------------


# process genotype data --------------------------------------------------------------------------
# genotypes here are a bit special here for F(n) intercross
# assuming one umi per genotype, each read tells us that underlying genotype is  either het or homozygous that genotype 
# assuming we had (more) depth data the emission probabilities defined here, informed by the 
# number of observed reads should factor into the posterior (https://github.com/tpbilton/GUSMap)
# what is the encoding here??
# dimensions are markers X cells
graw=readRDS(paste0(data.dir,'GenotypingPE6MiseqLibrary.rds'))

# approximate genetic map
# this is just a vector with same number of rows as graw 
gmap=readRDS(paste0(data.dir,'scRNAseqMiseqGeneticMap.rds'))
gmap.s=split_gmap(graw, gmap)

#if there is phenotype data, let's match everything up here! 
graw=graw[,match(rownames(dcounts), colnames(graw))]

#check if matrices match up by cell name
all.equal(colnames(graw), rownames(dcounts))
all.equal(colnames(graw), names(clusterID))
all.equal(colnames(graw), names(cell.counts))

cross=buildCrossObject(graw,gmap.s,'f2')

# run HMM (one iteration of forwards/backwards algorithm)
# get back probabilities for each of the different possible genotype combos
gps=calc_genoprob(cross, error_prob=.001)
#gps[[1]][1,,1]
#    AA     AB     BB 
#0.3563 0.5709 0.0728 

# additive allele dosage probabilities (effectively 0,1,2 mapped onto 0,.5,1)
gps.additive=genoprob_to_alleleprob(gps)
#A=AA+1/2(AB)
#B=BB+1/2(BB)
#R> gps.additive[[1]][1,,1]
#     A      B 
#0.6418 0.3582

# too many markers, just slowing things down, let's downsample to ~ 1 every cm
gps.s=reduce_markersJB(gmap.s,gps, 1)


#generate X matrix of additive covariates, annoying, let R do it for us
mm=model.matrix(dcounts[,1]~log2(cell.counts)+clusterID-1)
# what the heck does this look like?
# library(fields)
# image.plot(t(mm))

# hack for now
log2counts=log2(data.matrix(t(counts))+1)

# start with normal model and haley knott regression (here, includes dominance effect,
# can substitute with gps.additive to fit additive model, but looks like there's power
# loss)
gscan=scan1(gps.s,log2counts, addcovar=mm, cores=4)

# sanity check
# plot(which(gscan>5, arr.ind=T))
# lots of NAs coming from transcripts with 0 counts

# a few permutations------------------------
gscanPerm=list()
set.seed(10)
n.perms=5
for(i in 1:n.perms){
    print(i)
    snames=rownames(log2counts)
    s1=sample(1:nrow(log2counts))
    log2counts.p=log2counts[s1,]
    rownames(log2counts.p)=snames
    
    all.zeros=apply(log2counts.p, 2, function(x) sum(x==0)==nrow(log2counts))
    log2counts.p=log2counts.p[,-which(all.zeros)]

    addcovar.p=mm[s1,]
    rownames(addcovar.p)=snames
    gscanPerm[[i]]=scan1(gps.s,log2counts.p, addcovar=addcovar.p, cores=4)
}
#--------------------------------------------

# calculate FDR threshold, assumes exchangeability of max stat across transcripts
thresh=calcFDRthresh(gscan,gscanPerm)


# output table of peaks at FDR<10%
peaks=getPeaks(gscan,gmap.s,geneAnnot, thresh[2])

#plot peak positions
plot(peaks$peak.gcoord, peaks$t.gcoord)
abline(v=cumsum(sapply(split(geneAnnot$start_position,geneAnnot$chr), max)[-5]))
abline(h=cumsum(sapply(split(geneAnnot$start_position,geneAnnot$chr), max)[-5]))

#wtf, fix this 
ggplot(peaks, aes(x=peak_pmap_pos, y=start_position,size=lod))+
    geom_point()+
    ylab('transcript position')+xlab('QTL position')+
    facet_grid(rows=vars(chr),
               cols=vars(chromosome_name),
                scales='fixed', space='fixed',
                as.table=F, switch='both', shrink=T)    
