library(monocle3)
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

library(vcfR)
library(matrixStats)
library(spam)
# for parallelizing the HMM
ncores=24
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

sgd.genes=sgd.granges[sgd.granges$type=='gene',]
sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)


data.dir='/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/'
barcodes=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/barcodes.tsv'),sep='\t',header=F,stringsAsFactors=F)[,1]
counts=readMM(paste0(data.dir,'filtered_feature_bc_matrix/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]
features=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/features.tsv.gz'), sep='\t',header=F, stringsAsFactors=F)


aC=readMM(paste0(data.dir, 'alt_counts.mtx'))
rC=readMM(paste0(data.dir, 'ref_counts.mtx'))

variants=read.csv(paste0(data.dir,'out_var.txt'),sep='\t',header=F,stringsAsFactors=F)[,1]
variants= tstrsplit(variants,'_', type.convert=T)
variants[[2]]=variants[[2]]+1

all.vcf=read.vcfR(paste0(reference.dir, 'parents.nostar.vcf'))
gtsp=extract.gt(all.vcf)

#par.vcf=read.vcfR(paste0(reference.dir, 'BYxRMxYPS163xYJM145.vcf'))
par.vcf=read.vcfR(paste0(data.dir,'parents.vcf'))
gts=extract.gt(par.vcf)

#gts=gtsp[rownames(gts),c(2,9,11,16)]
gts[is.na(gts)]="0"

ee=.002

vname=paste0(variants[[1]],'_', variants[[2]])

rownames(rC)=vname
rownames(aC)=vname

#rownames(rC)=paste0(variants[[1]],'_', variants[[2]])
#colnames(rC)=barcodes

#rownames(aC)=rownames(rC)
#colnames(aC)=colnames(rC)
getLL=function(rC, aC, p1.ref, vname, ee=.002){
    #ivec=gts[,1]=='0'
    #p1.ref=gts[,1]=='0'

    p1=rbind(rC[p1.ref,], aC[!p1.ref,])
    #p1=as.dgCMatrix.spam(p1)
    vscramb=c(vname[p1.ref], vname[!p1.ref])
    rownames(p1)=vscramb
    p1=p1[vname,]
    p1=cleanup(as.spam.dgCMatrix(p1))

    p2=rbind(rC[!p1.ref,],aC[p1.ref,])
    #p2=as.dgCMatrix.spam(p2)
    vscramb=c(vname[!p1.ref], vname[p1.ref])
    rownames(p2)=vscramb
    p2=p2[vname,]
    p2=cleanup(as.spam.dgCMatrix(p2))

    n=p1+p2
    k=n-p2
   
    #rm(p1)
    #rm(p2)
    ps=dbinom(as.vector(n-k),as.vector(n),ee)
    #set ps==0 to 1 (bc log(1)=0), then can just do colsums
    ps[ps==0]=1
    lps=log(ps)

    pmat=n
    entries(pmat)=lps
    pmat=cleanup(pmat)

    ll=colSums(pmat)
    rm(ps)
    rm(n)
    rm(k)
    gc()
    
#    ll=apply.spam(pmat, 2, function(x) sum(log(x[x>0])))
    return(ll)
}  



#getLL=function(rC, aC, p1.ref,ee=.002){
#    #ivec=gts[,1]=='0'
#    #p1.ref=gts[,1]=='0'
#    p1=rbind(rC[p1.ref,],aC[!p1.ref,])
#    p1=p1[vname,]
#    
#    #norder=c(which(p1.ref),which(!p1.ref))
#    #p1=p1[norder,]
#
#    #p1=p1[rownames(rC),]
#    p2=rbind(rC[!p1.ref,],aC[p1.ref,])
#    #norder=c(which(!p1.ref),which(p1.ref))
#    #p2=p2[norder,]
#    p2=p2[vname,]
#    #p2=p2[rownames(rC),]
#
#    ref.counts=p1
#    alt.counts=p2
#
#    n=ref.counts+alt.counts
#    n=as.spam.dgCMatrix(n)
#    k=n-as.spam.dgCMatrix(alt.counts)
#    
#    ps=dbinom(as.vector(n-k),as.vector(n),ee)
#   
#    pmat=n
#    entries(pmat)=ps
#    pmat=cleanup(pmat)
#    ll=apply(pmat, 2, function(x) sum(log(x[x>0])))
#    return(ll)
#}  


BYLL  = getLL(rC,aC, gts[,1]=='0',vname)
RMLL  = getLL(rC,aC, gts[,2]=='0',vname)
YJMLL = getLL(rC,aC, gts[,3]=='0',vname)
YPSLL = getLL(rC,aC, gts[,4]=='0',vname)

cbind(BYLL, RMLL, YJMLL, YPSLL)->llik.table

m2diff=apply(llik.table,1, function(x) {
    y=sort(x,decreasing=T)
    (y[1]-y[2])
})

mcountl=apply(llik.table,1, which.max)
#mcountl[which(m2diff<585)]=NA
mcountl[mcountl==1]='BY'
mcountl[mcountl==2]='RM'
mcountl[mcountl==3]='YJM145'
mcountl[mcountl==4]='YPS163'

odir="/data/single_cell_eQTL/yeast/results/05_4-haploid-parental-strains/"
dir.create(odir)

saveRDS(mcountl, file=paste0(odir, "parental_assignment.RDS"))
saveRDS(llik.table, file=paste0(odir, "parental_llik.RDS"))

#mcountl2=readRDS(paste0(data.dir, "parental_assignment.RDS"))
library(monocle3)
mcountu=readRDS(paste0(odir, "parental_assignment.RDS"))
data.dir='/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/'
barcodes=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/barcodes.tsv'),sep='\t',header=F,stringsAsFactors=F)[,1]
counts=readMM(paste0(data.dir,'filtered_feature_bc_matrix/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]
features=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/features.tsv.gz'), sep='\t',header=F, stringsAsFactors=F)

genenames=features
names(genenames)[1]='name'
names(genenames)[2]='gene_short_name'

meta=data.frame(barcodes)
meta$parent=mcountl

cds=new_cell_data_set(counts,cell_metadata=meta, gene_metadata=genenames)
cds = preprocess_cds(cds, num_dim = 100)
#plot_pc_variance_explained(cds)
cds = reduce_dimension(cds)

plot_cells(cds, color_cells_by=as.factor(meta$parent), label_cell_groups=F)
u=(cds@reducedDims@listData$UMAP)
plot(u[,1], u[,2], pch=21, col='grey')
points(u[,1], u[,2],col=as.factor(meta$parent))
legend('topleft', legend=levels(as.factor(meta$parent)), fill=c(1,2,3,4))


sum(is.na(mcountl[which(u[,1]< (-6))]))
length(which(u[,1]< (-6)))


umap.file='/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/analysis/umap/2_components/projection.csv'
u=read_csv(umap.file)
plot(u[,1], u[,2], pch=21, col='grey')
points(u[,1], u[,2],col=as.factor(meta$parent))

