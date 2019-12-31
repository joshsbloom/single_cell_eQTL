library(monocle3)
library(Matrix)
library(rtracklayer)
library(data.table)
library(qtl)
library(qtl2)
library(parallel)
library(Rfast)
library(foreach)
library(doMC)
# for parallelizing the HMM
ncores=64
registerDoMC(cores=ncores)
#install.packages('doSNOW')
#library(doSNOW)
#cl=makeSOCKcluster(ncores)
#registerDoSNOW(cl)

# Data here from Eyal
# https://drive.google.com/drive/folders/12PlBkqHZgBZ_m8UK1upaL5EFDq5fJaPK?usp=sharing

# Data here from Josh (for running script and assessing output)
# https://drive.google.com/drive/folders/1fCr9Xc3NUin0RO6hnaEWXkIgRXouLD-g?usp=sharing

# various relevant data sets here
reference.dir='/data/single_cell_eQTL/elegans/reference/'

# additional functions here
code.dir='/data/single_cell_eQTL/elegans/code/'
source(paste0(code.dir, 'runHMM.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'utilities.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))


# where to output hmm results
hmm.out.dir='/data/single_cell_eQTL/elegans/XQTL_F4_2/hmm_v3/'
#vartrix count matrices at geno informative sites
vartrix.alt='/data/single_cell_eQTL/elegans/XQTL_F4_2/10xProcessed/altCounts_backgroundRemovedstrict.rds'
vartrix.ref='/data/single_cell_eQTL/elegans/XQTL_F4_2/10xProcessed/refCounts_backgroundRemovedstrict.rds'

# monocle data object 
monocle3.processed='/data/single_cell_eQTL/elegans/XQTL_F4_2/10xProcessed/backgroundCorrectedMonocle3_V2.rds'

monocle.data=readRDS(monocle3.processed)

#extract classification information-----------------
fine.classification=colData(monocle.data)@listData
counts=counts(monocle.data)
ncounts= normalized_counts(monocle.data)
#---------------------------------------------------

#add total counts to classification matrix----------
fine.classification$total=colSums(counts)
#---------------------------------------------------

#replace classifications with new object
#ctypes=readRDS('/data/single_cell_eQTL/elegans/reference/celltypeDataFrame_V3.rds')
#replace broad tissue classification 
#fine.classification$broad_tissue=ctypes@listData$broad_tissue

#new classifications
joint.covariates=data.frame(fine.classification)
#all.equal(names(fine.classification$Size_Factor), colnames(counts))
joint.covariates$barcode=colnames(counts)

all_broad=read.csv('/data/single_cell_eQTL/elegans/reference/classification/all_broad_pca_31oct.csv', stringsAsFactors=F)
doublets=read.csv('/data/single_cell_eQTL/elegans/reference/classification/doublets_all_31oct.csv', stringsAsFactors=F)
neuronal_pseudotimeE=read.csv('/data/single_cell_eQTL/elegans/reference/classification/neuronal_pseudotime_less_than50_31oct.csv', stringsAsFactors=F)
neuronal_pseudotime=read.csv('/data/single_cell_eQTL/elegans/reference/classification/neuronal_pseudotime_31oct.csv', stringsAsFactors=F)
joint.covariates$broad_tissue=NA
joint.covariates$fine_tissue=NA
joint.covariates$broad_tissue[match(all_broad$newid, joint.covariates$barcode)]=all_broad$broad_new
joint.covariates$PC1=NA
joint.covariates$PC2=NA
joint.covariates$PC1[match(all_broad$newid, joint.covariates$barcode)]=all_broad$PC1
joint.covariates$PC2[match(all_broad$newid, joint.covariates$barcode)]=all_broad$PC2

# using Eyal neuron classifications
joint.covariates$broad_tissue[match(neuronal_pseudotimeE$newid, joint.covariates$barcode)]='Neuron'
joint.covariates$fine_tissue[match(neuronal_pseudotimeE$newid, joint.covariates$barcode)]=neuronal_pseudotimeE$fine
joint.covariates$PC1[match(neuronal_pseudotimeE$newid, joint.covariates$barcode)]=neuronal_pseudotimeE$velocity_pseudotime


#get transcript annotations and put counts in genomic order -------------------
    transcript.data=rowData(monocle.data)
    transcript.annotations=import(paste0(reference.dir, 'genes.gtf'))
    cgenes=transcript.annotations[transcript.annotations$type=='gene',]
    # HACK!!, careful
    strand(cgenes)='*'
    cgenes=sort(cgenes)


    gcoord.key=c(0, cumsum(sapply(split(start(cgenes),seqnames(cgenes)), max)))[1:7]
    names(gcoord.key)[1:6]=chroms
    names(gcoord.key)[7]='MtDNA'

    transcript.data$gcoord=gcoord.key[transcript.data$chromosome_name]+transcript.data$start_position

    #rejigger counts matrces (order by position in genome)
    gene.reorder=na.omit(match(cgenes$gene_id, transcript.data$wbps_gene_id))
    #head(rownames(counts)[gene.reorder])

    counts=counts[gene.reorder,]
    ncounts=ncounts[gene.reorder,]
    transcript.data=transcript.data[gene.reorder,]
#------------------------------------------------------------------------------

# genetic map (this is for 10 generations) ------------------------------------
    gmap=readRDS(paste0(reference.dir, 'geneticMapXQTLsnplist.rds'))
    gmap$marker=paste0(gmap$chrom, '_', gmap$pos)
    #https://www.genetics.org/content/170/2/875
    #map size = X(AIL(J))=jX/2 (where j is number of generations); 
    #expected map size= 4*300/2 = 600cm; 1500cm/600cm = 2.5 (correction factor)
    # correction for F4
    gmap$map=gmap$map/2.5
#------------------------------------------------------------------------------

vartrix.parental.ref.file='/data/single_cell_eQTL/elegans/reference/refCountsL2.rds'
vartrix.parental.alt.file='/data/single_cell_eQTL/elegans/reference/altCountsL2.rds'
vartrix.parental.monocle.file='/data/single_cell_eQTL/elegans/reference/L2monocle3Object.rds'

# Run once -----------------------------------------------------------------------------
#this attaches N2.counts and CB.counts to the workspace ------------------------------
attach(getGenoInformativeCounts(vartrix.alt, vartrix.ref, gmap, 
                                      vartrix.parental.ref.file,
                                      vartrix.parental.alt.file,
                                      vartrix.parental.monocle.file  ))
#------------------------------------------------------------------------------

# extract info from genetic map from markers of interest
# jitter the map so transition probs from marker to marker are defined 
gmap.subset=gmap[match(rownames(N2.counts), gmap$marker),]
gmap.s =split(gmap.subset,gmap.subset$chrom)
#agh, don't want transition probs between markers to ever be 0, add some fudge factor between them
gmap.s =jitterGmap(gmap.s)

#calculate emission probabilities given counts
emissionProbs=estimateEmissionProbs(N2.counts, CB.counts)

## run the HMM
#runHMM(emissionProbs, hmm.out.dir, gmap.s, chroms,calc='genoprob')

## rerun and output viterbi path
#runHMM(emissionProbs, hmm.out.dir, gmap.s, chroms,calc='viterbi') 

## for historical/plotting/diagnosis purposes (hard calls of genotypes)
#rQTL.coded = encodeForRQTL_F2(N2.counts, CB.counts)

##previous r/qtl HMM ouput
#G=readRDS('/data/single_cell_eQTL/elegans/XQTL_F4_2/additive.geno_v2.RDS')
#G=G[,colnames(G) %in%  rownames(N2.counts)]

# load everything back in, and sanity check
posteriorProbs=list()
viterbiPath=list()
additiveCoding=list()
for(cc in chroms){
    print(cc)
    posteriorProbs[[cc]]=readRDS(paste0(hmm.out.dir, cc,'.RDS'))
    # viterbiPath[[cc]]=readRDS(paste0(hmm.out.dir, cc,'_viterbi.RDS'))
    #probability of CB allele
    additiveCoding[[cc]]=.5*posteriorProbs[[cc]][,2,]+posteriorProbs[[cc]][,3,]
}
geno=do.call('cbind', additiveCoding)
saveRDS(geno, file='/data/single_cell_eQTL/elegans/results/20191217/additiveGenoCoding.RDS')
# -------------------------------------------------------------------------------Run Once -----
##G=do.call('cbind', additiveCoding)

#for( k in 13000:15000){
#    diagnoseHMM(k, chrom='II', gmap.s, N2.counts,CB.counts, additiveCoding, rQTL.coded)#viterbiPath=NULL,  G=NULL)
#    readline()
#}

#save this for speedup
#v=do.call('cbind', viterbiPath)

rm(posteriorProbs)
rm(additiveCoding)
rm(viterbiPath)
gc()
