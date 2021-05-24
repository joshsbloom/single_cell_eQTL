library(Matrix)
library(data.table)
library(qtl)
library(qtl2)
library(parallel)
library(Rfast)
library(seqinr)
library(GenomicRanges)
#library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(foreach)
library(doMC)
library(glmmTMB)
library(ggplot2)
library(vcfR)
library(monocle)
library(matrixStats)
library(spam)


library(foreach)
library(doMC)
library(dplyr)
# for parallelizing the HMM

code.dir='/data/single_cell_eQTL/yeast/code/'
reference.dir='/data/single_cell_eQTL/yeast/reference/'
data.base.dir='/data/single_cell_eQTL/yeast/processed/'

source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))

#load custom functions for processing ASE data 
source(paste0(code.dir, 'ASE_fxs.R'))

sgd.genes=getSGD_Gene_Intervals(reference.dir)

experiment.name='12_Group_1_diploids_3004_2444_3051_5k_Feb_21/'

#order matters here for counting diploid specific variants, so proceed carefully
input.diploids=names(crosses.to.parents)[c(2,4,14)]

data.dir=paste0(data.base.dir, experiment.dir)

#vname=paste0(variants[[1]],'_', variants[[2]])
ase.Data=buildASE_data(data.dir)

dip.Assignments=getDipAssignments(ase.Data, input.diploids, crosses.to.parents, ncores=8)
dip.specificCounts=countDiploidSpecificVariants(ase.Data, input.diploids, crosses.to.parents)
#overwrite dip.Assignments to include likelihood, counts, and umap
dip.Assignments=dplyr::left_join(dplyr::left_join(ase.Data$um, dip.Assignments, by='barcode'), dip.specificCounts)

for(dip in input.diploids) {
    dip='A'
    phasedCounts=getPhasedCountsPerTranscript(ase.Data, dip.Assignments, dip, sgd.genes, crosses.to.parents) 
}







rMat=phasedCounts[[1]] #ref.ASEcounts #do.call('cbind', ref.ASEcounts)
aMat=phasedCounts[[2]] #alt.ASEcounts # do.call('cbind', alt.ASEcounts)

plot(log2(rowSums(rMat)), log2(rowSums(aMat)))
abline(0,1)

plot(log2(colSums(rMat)), log2(colSums(aMat)))
abline(0,1)

totRcounts=rowSums(rMat)
totAcounts=rowSums(aMat)

totGcounts=colSums(rMat+aMat)

A_R=colSums(aMat)/colSums(rMat)
aseElife=gdata::read.xls('/data/single_cell_eQTL/yeast/reference/elife-35471-data7-v2_ASE_results.xlsx')
raw.ratio=data.frame(gene=names(A_R), ratio=as.vector(A_R),sum=totGcounts)
#library(tidyverse)
test=left_join(aseElife, raw.ratio, by='gene')
test=test[test$ratio!=0 & is.finite(test$ratio),]

plot(log2(test$ratio),test$log2ASEFoldChange)



pairs(pae[,12:14],col=pae$diploid_assignment)
pairs(pae[,5:7],col=pae$diploid_assignment)
pairs(pae[,2:3],col=pae$diploid_assignment)

#diagnostic plots 
par(mfrow=c(1,3))
with(pae,{
plot(rcount.A, rcount.B, col=diploid_assignment,
     xlab='by x rm rare variant counts', ylab='ypsxyjm rare variant counts')
plot(rcount.A,rcount.3004, col=diploid_assignment,
     xlab='by x rm rare variant counts', ylab='3004 rare variant counts')
plot(rcount.B,rcount.3004, col=diploid_assignment,
     xlab='yps x yjm rare variant counts', ylab='3004 rare variant counts'
)
})


library(ggplot)
library(viridis)
library(ggpubr)
a=ggplot(pae, aes(x=UMAP.1, y=UMAP.2, color=log10(umis_per_cell))) + 
    geom_point() + scale_color_viridis()
b=ggplot(pae, aes(x=UMAP.1, y=UMAP.2, color=diploid) ) +
   geom_point() 
#cc=ggplot(um, aes(x=UMAP.1, y=UMAP.2, color=diploid_assignment_diff) ) +
#   geom_point() + scale_color_viridis()

cc=ggplot(pae, aes(x=UMAP.1, y=UMAP.2, color=(rcount.3004)) ) +
   geom_point() + scale_color_viridis()+
   ggtitle('colored by YJM981 x CBS2888a specific variant counts per cell')
dd=ggplot(pae, aes(x=UMAP.1, y=UMAP.2, color=(rcount.B)) ) +
   geom_point() + scale_color_viridis()+
   ggtitle('colored by YJM145x x YPS163a specific variant counts per cell')
ee=ggplot(pae, aes(x=UMAP.1, y=UMAP.2, color=(rcount.A)) ) +
   geom_point() + scale_color_viridis()+
   ggtitle('colored by BY x RM specific variant counts per cell')

ggarrange(a,b, cc, dd, ee, nrow=2, ncol=3)






























#nn=colnames(gts)[1]
#p1.ref=gts[,nn]=='0'

#lls=list()
# be careful about memory usage here 


#mcountl[m2diff<500]=NA

#r-evaluate filters (for fitst diploid panel it was 585)
#mcountl[m2diff<585]=NA
#mcountl[m2diff<400]=NA
saveRDS(mcountl, file=paste0(data.dir, "parental_assignment.RDS"))

#counts at geno informative sites 
#gcount=(rC+aC)

countM=readMM(paste0(data.dir,'filtered_feature_bc_matrix/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]
tcount=colSums(countM)
acount=apply(gts, 1, function(x) sum(x=='1'))
rvars=acount==1
gts.sub=gts[rvars,]
ap1=(gts.sub[,1]=="1" | gts.sub[,2]=="1")
ap2=(gts.sub[,3]=="1" | gts.sub[,4]=="1")
ap3=(gts.sub[,5]=="1" | gts.sub[,6]=="1")

aCr=aC[rvars,]
ap1c=colSums(aCr[which(ap1),]>0)
ap2c=colSums(aCr[which(ap2),]>0)
ap3c=colSums(aCr[which(ap3),]>0)



um=read_csv(paste0(data.dir, 'analysis/umap/2_components/projection.csv'))
plot(um[,1], um[,2], col=mcountl)



um$umis_per_cell=colSums(countM)
um$diploid_assignment=mcountl
um$diploid = as.vector(sapply(crosses.to.parents[colnames(dipLik)], paste, collapse=' x '))[mcountl]
um$diploid_assignment_likdiff=m2diff
um$ap1c=ap1c
um$ap2c=ap2c
um$ap3c=ap3c



um[,1]=as.numeric(um[,1])
um[,2]=as.numeric(um[,2])
names(um)[1:2]=c('umap1', 'umap2')
names(um)[7:9]=colnames(dipLik)

saveRDS(um, file=paste0(data.dir, 'parental_assignment_extended.RDS'))





#uBYxRM=gts[,2]!=gts[,1] & gts[,2]!=gts[,3] & gts[,2]!=gts[,4] & gts[,2] != gts[,5] & gts[,2]!=gts[,6]
#uYPSxYJM.1=(gts[,3]!=gts[,1]) & gts[,3]!=gts[,2] & gts[,3]!=gts[,4] & gts[,3] != gts[,5] & gts[,3]!=gts[,6]
#uYPSxYJM.2=(gts[,4]!=gts[,1]) & gts[,4]!=gts[,2] & gts[,4]!=gts[,3] & gts[,4] != gts[,5] & gts[,4]!=gts[,6]
#u3004=
#aC[which(gts[,2]!=gts[,1] & gts[,2]!=gts[,3] & gts[,2]!=gts[,4] & gts[,2] != gts[,5] & gts[,2]!=gts[,6]),]



#rownames(rC)=paste0(variants[[1]],'_', variants[[2]])
#colnames(rC)=barcodes

#rownames(aC)=rownames(rC)
#colnames(aC)=colnames(rC)
#getLL=function(rC, aC, p1.ref,ee=.002){
#    
#    p1=rbind(rC[p1.ref,],aC[!p1.ref,])
#    norder=c(which(p1.ref),which(!p1.ref))
#    p1=p1[norder,]
#
#    #p1=p1[rownames(rC),]
#    p2=rbind(rC[!p1.ref,],aC[p1.ref,])
#    norder=c(which(!p1.ref),which(p1.ref))
#    p2=p2[norder,]
#    #p2=p2[rownames(rC),]
#
#    #total
#    n=p1+p2
#
#    # just p1 
#    k=n-p2
#    rm(p1); rm(p2);
#    gc()
#    # compute likelihoods    
#    ps=dbinom(as.vector(n-k),as.vector(n),ee,log=T)
#    pmat=n
#    entries(pmat)=ps
#    pmat=cleanup(pmat)
#    #return(pmat)
#    rm(n); rm(k);
#    ll=colSums(pmat)
#    return(ll)
#}  

















tm=table(mcountl)
names(tm)=colnames(dipLik)
barplot(tm)

counts=readMM(paste0(data.dir,'filtered_feature_bc_matrix/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]

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






#some genome annotation
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

sgd.genes=sgd.granges[sgd.granges$type=='gene',]
sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)


#data.dir='/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/'
data.dir='/data/single_cell_eQTL/yeast/processed/13_Group1_diploids_3004_2444_3051_7point5K_Feb_21/'

#data.dir='/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/'
#rownames(rC)=vname
#rownames(aC)=vname


