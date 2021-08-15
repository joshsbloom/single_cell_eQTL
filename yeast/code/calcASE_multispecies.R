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
library(monocle)
library(matrixStats)
library(spam)
library(dplyr)
library(broom.mixed)

code.dir='/data/single_cell_eQTL/yeast/code/'
source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))


chroms=paste0('chr', as.roman(1:16)) 
#c('I','II', 'III', 'IV', 'V', 'X')

#some genome annotation --------------------------------------------------------
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
#sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
#sgd=as.data.frame(sgd.granges)
#gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

#sgd.genes=sgd.granges[sgd.granges$type=='gene',]
#sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)
#--------------------------------------------------------------------------------

#update gene definitions with utr
sgd.granges=import.gff(paste0(reference.dir, 'genes.gtf'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))
sgd.genes=sgd.granges[sgd.granges$type=='exon',]
sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)



#data.dir='/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/'
#data.dir='/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/'
#data.dir='/data/single_cell_eQTL/yeast/processed/12_Group_1_diploids_3004_2444_3051_5k_Feb_21/'

#multi-species data here : 
#https://drive.google.com/drive/u/1/folders/1eWEzKtzA2sddibi9tOpuxPZi7an7e4EM

best_orthologs =read.csv('/data/single_cell_eQTL/yeast/results/multi_species/best_orthologs.tsv',sep='\t',header=T,stringsAsFactors=F)


dip.oi='SCER_SEUB'
dip.oi='SPAR_SEUB'

assignment.e1 = read.csv("/data/single_cell_eQTL/yeast/results/multi_species/GP-2-5k_SCER_and_SEUB_and_SPAR/cell_assignment.txt",sep='\t',header=T,stringsAsFactors=F)
data.dir.e1='/data/single_cell_eQTL/yeast/results/multi_species/GP-2-5k_SCER_and_SEUB_and_SPAR/'
assignment.e2 = read.csv("/data/single_cell_eQTL/yeast/results/multi_species/GP-2-7point5k_SCER_and_SEUB_and_SPAR/cell_assignment.txt",sep='\t',header=T,stringsAsFactors=F)
data.dir.e2='/data/single_cell_eQTL/yeast/results/multi_species/GP-2-7point5k_SCER_and_SEUB_and_SPAR/'

cc.e1=read_csv(paste0('/data/single_cell_eQTL/yeast/results/multi_species/GP-2-5k_SCER_and_SEUB_and_SPAR/', dip.oi, '/cell_cycle_assignments.csv'))
cc.e2=read_csv(paste0('/data/single_cell_eQTL/yeast/results/multi_species/GP-2-7point5k_SCER_and_SEUB_and_SPAR/', dip.oi, '/cell_cycle_assignments.csv'))


#data.dir='/data/single_cell_eQTL/yeast/processed/18_3051_May10/'
#2_Group_1_diploids_3004_2444_3051_5k_Feb_21/'

barcodes.e1=read.csv(paste0(data.dir.e1,'barcodes.tsv'),
                  sep='\t',header=F,stringsAsFactors=F)[,1]
barcodes.e2=read.csv(paste0(data.dir.e2,'barcodes.tsv'),
                  sep='\t',header=F,stringsAsFactors=F)[,1]

features.e1=read.csv(paste0(data.dir.e1,'features.tsv'), 
                  sep='\t',header=F, stringsAsFactors=F)
features.e2=read.csv(paste0(data.dir.e2,'features.tsv'), 
                  sep='\t',header=F, stringsAsFactors=F)

counts.e1=spam::read.MM(paste0(data.dir.e1,'/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]
counts.e2=spam::read.MM(paste0(data.dir.e2,'/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]

numi.e1=colSums(counts.e1)
numi.e2=colSums(counts.e2)

#data.dir='/data/single_cell_eQTL/yeast/processed/18_3051_May10/'
#2_Group_1_diploids_3004_2444_3051_5k_Feb_21/'

to.keep.e1=assignment.e1$cell_id[assignment.e1$best_match==dip.oi & assignment.e1$pass=="Pass" & assignment.e1$cell_id %in% cc.e1$cell_name]
to.keep.e2=assignment.e2$cell_id[assignment.e2$best_match==dip.oi & assignment.e2$pass=="Pass" & assignment.e2$cell_id %in% cc.e2$cell_name]

cc.e1=cc.e1[cc.e1$cell_name %in% to.keep.e1,]
cc.e2=cc.e2[cc.e2$cell_name %in% to.keep.e2,]

#all.equal(barcodes.e1[barcodes.e1 %in% to.keep.e1], cc.e1$cell_name)

numi.e1=numi.e1[which(barcodes.e1 %in% to.keep.e1)]
numi.e2=numi.e2[ which(barcodes.e2 %in% to.keep.e2)]

#bsub.e1=barcodes.e1[barcodes.e1 %in% to.keep.e1]

# for scer replace first - only
if(dip.oi=='SCER_SEUB'){
csub.e1.1=counts.e1[match(sub('-', '_', best_orthologs$gene_id.scer), features.e1[,2]), which(barcodes.e1 %in% to.keep.e1)]
fsub.e1.1=features.e1[match(sub('-', '_', best_orthologs$gene_id.scer), features.e1[,2]),]
csub.e2.1=counts.e2[match(sub('-', '_', best_orthologs$gene_id.scer), features.e2[,2]), which(barcodes.e2 %in% to.keep.e2)]
fsub.e2.1=features.e2[match(sub('-', '_', best_orthologs$gene_id.scer), features.e2[,2]),]
}
if(dip.oi=='SPAR_SEUB'){
    csub.e1.1=counts.e1[match(gsub('-', '_', best_orthologs$gene_id.spar), features.e1[,2]), which(barcodes.e1 %in% to.keep.e1)]
    fsub.e1.1=features.e1[match(gsub('-', '_', best_orthologs$gene_id.spar), features.e1[,2]),]
    csub.e2.1=counts.e2[match(gsub('-', '_', best_orthologs$gene_id.spar), features.e2[,2]), which(barcodes.e2 %in% to.keep.e2)]
    fsub.e2.1=features.e2[match(gsub('-', '_', best_orthologs$gene_id.spar), features.e2[,2]),]
}




# for seub replace all - 
csub.e1.2=counts.e1[match(gsub('-', '_', best_orthologs$gene_id.seub), features.e1[,2]),which(barcodes.e1 %in% to.keep.e1)]
fsub.e1.2=features.e1[match(gsub('-', '_', best_orthologs$gene_id.seub), features.e1[,2]),]
csub.e2.2=counts.e2[match(gsub('-', '_', best_orthologs$gene_id.seub), features.e2[,2]),which(barcodes.e2 %in% to.keep.e2)]
fsub.e2.2=features.e2[match(gsub('-', '_', best_orthologs$gene_id.seub), features.e2[,2]),]


non.zero.cells.e1=apply((csub.e1.1 +csub.e1.2), 1, function(x) sum(x>0))
non.zero.cells.e2=apply((csub.e2.1 +csub.e2.2), 1, function(x) sum(x>0))

genes.to.keep=(non.zero.cells.e1+non.zero.cells.e2)>128 #  & non.zero.cells.e2>64 # numi[paroi]<10000

csub.e1.1=csub.e1.1[genes.to.keep,]
csub.e1.2=csub.e1.2[genes.to.keep,]

csub.e2.1=csub.e2.1[genes.to.keep,]
csub.e2.2=csub.e2.2[genes.to.keep,]


non.zero.cells.p1=apply(csub.e1.1, 1, function(x) sum(x>0))+ apply(csub.e2.1, 1, function(x) sum(x>0)) 
non.zero.cells.p2=apply(csub.e1.2, 1, function(x) sum(x>0))+ apply(csub.e2.2, 1, function(x) sum(x>0)) 
tot.counts.p1=apply(csub.e1.1, 1, function(x) sum(x))+ apply(csub.e2.1, 1, function(x) sum(x)) 
tot.counts.p2=apply(csub.e1.2, 1, function(x) sum(x))+ apply(csub.e2.2, 1, function(x) sum(x)) 

stats=data.frame(gene_id.scer=best_orthologs$gene_id.scer[genes.to.keep],
           non.zero.cells.p1=non.zero.cells.p1,
           non.zero.cells.p2= non.zero.cells.p2,
           tot.counts.p1=tot.counts.p1,
           tot.counts.p2= tot.counts.p2,
           tot.cells=ncol(csub.e1.1)+ncol(csub.e2.1)
            )
saveRDS(stats, file=paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/stats_', dip.oi,'.RDS'))
           

#fsub.1=fsub.1[genes.to.keep,]
#fsub.2=fsub.2[genes.to.keep,]

best_orthologs.to.keep=best_orthologs[genes.to.keep,]

#beta binomial test for ASE
bbin=mclapply(1:nrow(csub.e1.1),
    function(g,...){
        r=c(as.vector(as.matrix(csub.e1.1[g,])),
            as.vector(as.matrix(csub.e2.1[g,])))
        a=c(as.vector(as.matrix(csub.e1.2[g,])),
            as.vector(as.matrix(csub.e2.2[g,])))
        #efs=experiment.factor[paroi][numi[paroi]<5000]
        ex=cbind(a,r)#[numi[paroi]<5000,] #a)
        #dcount[[g]]=colSums(ex)
        ef=c(rep('E1', ncol(csub.e1.1)),rep('E2', ncol(csub.e2.1)))
        print(g)
        return(test=glmmTMB(ex~ef,family=glmmTMB::betabinomial(link='logit'))) 
    },
    csub.e1.1,csub.e1.2,csub.e2.1,csub.e2.2,mc.cores=48)


names(bbin)=best_orthologs$gene_id.scer[genes.to.keep]
bbin.llik=sapply(bbin, logLik)
bbin.model.converged.genes=names(bbin.llik)[!is.na(bbin.llik)]

#why so slow??????
lbb=lapply(bbin[bbin.model.converged.genes], function (x) coef(summary(x))$cond)
rm(bbin)
lbbdf=rbindlist(lapply(lbb, function(x) data.frame(x)[1,]), idcol='gene')
saveRDS(lbbdf, file=paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/bbin_', dip.oi, '.RDS'))

#write_csv(scASE_18, file='~/Dropbox/FileTransfer/scASE_18.csv')
#'/home/jbloom/Dropbox/FileTransfer/single_cell/multispecies/bbin_SCERSEUB_llikCC.RDS')




bbinCC=mclapply(1:nrow(csub.e1.1),
    function(g,...){
        r=c(as.vector(as.matrix(csub.e1.1[g,])),
            as.vector(as.matrix(csub.e2.1[g,])))
        a=c(as.vector(as.matrix(csub.e1.2[g,])),
            as.vector(as.matrix(csub.e2.2[g,])))
        cc=c(as.vector(cc.e1$cell_cycle), as.vector(cc.e2$cell_cycle))
        #efs=experiment.factor[paroi][numi[paroi]<5000]
        ex=cbind(a,r)#[numi[paroi]<5000,] #a)
        #dcount[[g]]=colSums(ex)
        ef=c(rep('E1', ncol(csub.e1.1)),rep('E2', ncol(csub.e2.1)))
        print(g)
        cc=as.factor(cc)
        cc=relevel(cc, 'M/G1')
        return(test=glmmTMB(ex~ef+cc,family=glmmTMB::betabinomial(link='logit'))) 
    },
    csub.e1.1,csub.e1.2,csub.e2.1,csub.e2.2, cc.e1, cc.e2, mc.cores=48)
names(bbinCC)=best_orthologs$gene_id.scer[genes.to.keep]
bbin.llikCC=sapply(bbinCC, logLik)
bbin.model.converged.genesCC=names(bbin.llikCC)[!is.na(bbin.llikCC)]


lbbCC=lapply(bbinCC[bbin.model.converged.genesCC], function (x) coef(summary(x))$cond)
rm(bbinCC)
saveRDS(lbbCC, file=paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/bbin_', dip.oi, '_CC.RDS'))
lbbdfCC=rbindlist(lapply(lbbCC, function(x) data.frame(x)[1,]), idcol='gene')
binCC.meltl=lapply(lbbCC, function(x) melt(x) )
binCC.meltl=rbindlist(binCC.meltl, idcol='gene')
saveRDS(binCC.meltl, paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/bbin_', dip.oi, 'CC_melted.RDS'))


llik.bbin.cc=data.frame(bbin.llik=bbin.llik, bbbinCC.llik=bbin.llikCC, LRS=-2*(bbin.llik-bbin.llikCC), p.value=pchisq(-2*(bbin.llik-bbin.llikCC), 4,lower.tail=F))
saveRDS(llik.bbin.cc, paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/bbin_', dip.oi, 'llikCC.RDS') )





scASE_18=lbbdf

ncells=sum(c(ncol(csub.e1.1),ncol(csub.e2.1)))
cellID=factor(as.character(seq(1:ncells))) #length(r)

setup.vars=list(
        cellIDf=c(cellID,cellID),
        geno=factor(c(rep('A',ncells),rep('B',ncells))),
        of2=c(log(numi.e1),log(numi.e2), log(numi.e1), log(numi.e2)),
        ef=c(rep('E1', ncol(csub.e1.1)),rep('E2', ncol(csub.e2.1)),
             rep('E1', ncol(csub.e1.1)),rep('E2', ncol(csub.e2.1)))
       )

#cl <- makeCluster(48)
#clusterEvalQ(cl, {
#       library(glmmTMB)
#       #library(Matrix)
#       #library(MASS) 
#       NULL  })
#clusterExport(cl, varlist=c("thetas", "Yks", "Gsub", "mmp1", "nbLL", "domap"))
#clusterEvalQ(cl,{ Y=Yks;    DM=mmp1;   return(NULL);})


nbin=list()
nbin=mclapply(1:nrow(csub.e1.1), 
   function(g, ... ){
        print(g)
        r=c(as.vector(as.matrix(csub.e1.1[g,])),
            as.vector(as.matrix(csub.e2.1[g,])))
        a=c(as.vector(as.matrix(csub.e1.2[g,])),
            as.vector(as.matrix(csub.e2.2[g,])))

      #  r=rMat[cells.to.keep,g]
      #  a=aMat[cells.to.keep,g]
        #efs=experiment.factor[paroi][numi[paroi]<5000]
        y=c(r,a)#[numi
        cellIDf=setup.vars$cellIDf
        of2=setup.vars$of2
        geno=setup.vars$geno
        ef=setup.vars$ef
       # (1|cellIDf)
        nbfit=glmmTMB(y~geno+ef+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
        if(!is.na(logLik(nbfit))){     nbin[[g]]=summary(nbfit)  } else {nbin[[g]]=NULL }
   },
    csub.e1.1,csub.e1.2,csub.e2.1,csub.e2.2,setup.vars,mc.cores=48)
# rMat,aMat,cells.to.keep,setup.vars,mc.cores=48)
names(nbin)=best_orthologs$gene_id.scer[genes.to.keep]

nbin.model.converged.genes=names(nbin)[sapply(nbin, function(x) !is.null(x))]
nbin.reduced=lapply(nbin[nbin.model.converged.genes], function(x) list(coef=x$coef, stats=x$AICtab))
saveRDS(nbin.reduced, file=paste0('/data/single_cell_eQTL/yeast/results/multi_species/results/nbin_', dip.oi, 'raw.RDS'))
































#sapply(nbin.reduced, function(x) x$stats)
#
#nbin.disp=list()
#nbin.disp=mclapply(1:nrow(csub.e1.1),  
#   function(g, ... ){
#        print(g)
#        r=c(as.vector(as.matrix(csub.e1.1[g,])),
#            as.vector(as.matrix(csub.e2.1[g,])))
#        a=c(as.vector(as.matrix(csub.e1.2[g,])),
#            as.vector(as.matrix(csub.e2.2[g,])))
#
#      #  r=rMat[cells.to.keep,g]
#      #  a=aMat[cells.to.keep,g]
#        #efs=experiment.factor[paroi][numi[paroi]<5000]
#        y=c(r,a)#[numi
#        cellIDf=setup.vars$cellIDf
#        of2=setup.vars$of2
#        geno=setup.vars$geno
#        ef=setup.vars$ef
#      #  nbfitg=glmmTMB(y~geno+ef+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
#        nbfit=glmmTMB(y~ef+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
#        if(!is.na(logLik(nbfit))){     nbin.disp[[g]]=summary(nbfit)  } else {nbin.disp[[g]]=NULL }
#   },
#    csub.e1.1,csub.e1.2,csub.e2.1,csub.e2.2,setup.vars,mc.cores=48)
#names(nbin.disp)=best_orthologs$gene_id.scer[genes.to.keep]
#
#nbin.disp.model.converged.genes=names(nbin.disp)[sapply(nbin.disp, function(x) !is.null(x))]
#nbin.disp.reduced=lapply(nbin.disp[nbin.disp.model.converged.genes], function(x) list(coef=x$coef, stats=x$AICtab))
#
#
#saveRDS(nbin.disp.reduced, file='/data/single_cell_eQTL/yeast/results/multi_species/results/nbin_nogeno_SCERSEUBraw.RDS')
#
#f=sapply(nbin.reduced, function(x) x$stats)
#r=sapply(nbin.disp.reduced, function(x) x$stats)
#
#f=data.frame(t(f),gene=colnames(f))
#r=data.frame(t(r),gene=colnames(r))
#
#fr=(left_join(f,r, by='gene'))
#
#
#
#nbin.model.converged.genes=names(nbin)[sapply(nbin, function(x) !is.null(x))]
#scASE_18_nbin.mean=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$cond[2,])))
#scASE_18_nbin.varInt=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$disp[1,])))
#scASE_18_nbin.var=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$disp[2,])))
#
#scASE_18_nbin=left_join(scASE_18_nbin.mean, scASE_18_nbin.var, by="gene",suffix=c(".mean", ".var"))
#saveRDS(scASE_18_nbin, file='/data/single_cell_eQTL/yeast/results/multi_species/results/nbin_SCERSEUB.RDS')
#
#
#
#mElife_scASE_18=left_join(aseElife, scASE_18, by='gene')
#bbin_nbin=left_join(scASE_18, scASE_18_nbin,by='gene')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#variants=read.csv(paste0(data.dir,'out_var.txt'),sep='\t',header=F,stringsAsFactors=F)[,1]
#variants= tstrsplit(variants,'_', type.convert=T)
#variants[[2]]=variants[[2]]+1
#
##cc=read.csv(paste0(data.dir,'cell_cycle_3051_4.txt'),sep=' ',header=T,stringsAsFactors=F)
#
#markerGR=makeGRangesFromDataFrame(data.frame(chr=variants[[1]], 
#                                             start=as.numeric(variants[[2]]),
#                                             end=as.numeric(variants[[2]]), 
#                                             strand="*"))
#markerGR$name=paste0(variants[[1]], ':', variants[[2]]) 
#markerGR$gcoord=gcoord.key[as.character(seqnames(markerGR))]+start(markerGR)
#markerGR$name2=paste0(variants[[1]], '_', variants[[2]])
##rownames(aC)=markerGR$name2
##rownames(rC)=markerGR$name2
#
#gts.subset=gts[rownames(gts)%in%markerGR$name2,]
#
#
#bad.markers=which(!markerGR$name2%in%rownames(gts.subset))
##careful, assuming bad.markers is not empty
#markerGR=markerGR[-bad.markers,]
#aC=aC[-bad.markers,]
#rC=rC[-bad.markers,]
#
#ee=.002
##p.assign=readRDS(paste0(data.dir, "parental_assignment.RDS"))
#
#p.assign=readRDS(paste0(data.dir, "parental_assignment_extended.RDS")) #$diploid_assignment
##hack
##p.assign=rep(2, ncol(aC))
#
#to.keep=which(!is.na(p.assign))
#numiG=colSums(aC)+colSums(rC)
#
#cross.index=2
#paroi=which(p.assign$diploid=="BYa x RMx" ) #cross.index)
##paroi=which(p.assign$diploid=="YJM981x x CBS2888a" ) #cross.index)
#
#rCs=rC[,paroi]
#aCs=aC[,paroi]
#
## rejigger logic , 5/10/21
#cross.index=2
##cross.index=14
##crosses.to.parents[[cross.index]]
#p1.ref=gts.subset[, crosses.to.parents[[cross.index]][1]]=='0' 
#p1.not.p2=gts.subset[, crosses.to.parents[[cross.index]][2]]!=gts.subset[, crosses.to.parents[[cross.index]][1]]
#
##p1.refwhich(is.na(p1.ref))]=TRUE
#
##pull out counts at all the sites corresponding to parent 1 
#p1=rbind(rCs[p1.ref & p1.not.p2,],aCs[!p1.ref & p1.not.p2 ,])
#p1=as.matrix(p1)
#rownames(p1)=names(c(which(p1.ref & p1.not.p2),which(!p1.ref & p1.not.p2)))
##norder=c(which(p1.ref),which(!p1.ref))
##p1=p1[norder,]
#
##pull out counts at all the sites corresponding to parent 2
#p2=rbind(rCs[!p1.ref & p1.not.p2,],aCs[p1.ref & p1.not.p2,])
#p2=as.matrix(p2)
#rownames(p2)=names(c(which(!p1.ref& p1.not.p2),which(p1.ref & p1.not.p2)))
#
##norder=c(which(!p1.ref),which(p1.ref))
##p2=p2[norder,]
#
##rownames(p1)=markerGR$name
#    #p2=p2[rownames(rC),]
#g.counts=list(p1,p2)
#
#m.granges=markerGR[markerGR$name2 %in% rownames(p1)]
#
#
#fO=findOverlaps(sgd.genes, m.granges)
##mgsplit=split(m.granges$name[subjectHits(fO)], sgd.genes$Name[queryHits(fO)])
#mgsplit=split(m.granges$name2[subjectHits(fO)], sgd.genes$gene_id[queryHits(fO)])
#
#print('calc counts per gene, allele 1')
#
## could restructure these as sparse matrices if necessary
#threads=32
#ref.ASEcounts=mcmapply( function(x) {
#            #m=match(x, markerGR$name2)                               
#            if(length(x)>1) {
#                colSums(g.counts[[1]][x,])
#            } else {
#                g.counts[[1]][x,]
#            }
#       },mgsplit,mc.cores=threads)
#
#print('calc counts per gene, allele 2')
#alt.ASEcounts=mcmapply( function(x) {
#            #m=match(x, markerGR$name2)                               
#            if(length(x)>1) {
#                colSums(g.counts[[2]][x,])
#            } else {
#                g.counts[[2]][x,]
#            }
#}, mgsplit,mc.cores=threads)
#
#rm(g.counts)
#rm(p1)
#rm(p2)
#
##r.lengths=sapply(ref.ASEcounts, length)
##a.lengths=sapply(alt.ASEcounts, length)
#
##weird way of pulling out the genes that don't have counts in all the cells
##bgenes=unique(c( names(which(a.lengths< max(r.lengths))), names(which(r.lengths< max(r.lengths)))))
##ref.ASEcounts=ref.ASEcounts[!names(ref.ASEcounts)%in%bgenes]
##alt.ASEcounts=alt.ASEcounts[!names(alt.ASEcounts)%in%bgenes]
#rMat=ref.ASEcounts #do.call('cbind', ref.ASEcounts)
#aMat=alt.ASEcounts # do.call('cbind', alt.ASEcounts)
#
#plot(log2(rowSums(rMat)), log2(rowSums(aMat)))
#abline(0,1)
#
#plot(log2(colSums(rMat)), log2(colSums(aMat)))
#abline(0,1)
#
#totRcounts=rowSums(rMat)
#totAcounts=rowSums(aMat)
#
#totGcounts=colSums(rMat+aMat)
#
#A_R=colSums(aMat)/colSums(rMat)
#aseElife=gdata::read.xls('/data/single_cell_eQTL/yeast/reference/elife-35471-data7-v2_ASE_results.xlsx')
#raw.ratio=data.frame(gene=names(A_R), ratio=as.vector(A_R),sum=totGcounts)
##library(tidyverse)
#test=left_join(aseElife, raw.ratio, by='gene')
#test=test[test$ratio!=0 & is.finite(test$ratio),]
#
#plot(log2(test$ratio),test$log2ASEFoldChange)
#cor.test(test$log2ASEFoldChange, log2(test$ratio), method='spearman')
#
#plot(test$log2ASEFoldChange[log2(test$sum)>10], log2(test$ratio)[log2(test$sum)>10])
#cor.test(test$log2ASEFoldChange[log2(test$sum)>5], log2(test$ratio)[log2(test$sum)>5])
#
#plot(log2(rowSums(rMat)[numi[paroi]<5000]), log2(rowSums(aMat)[numi[paroi]<5000]))
#
#
##some thresholding on what a cell is, if more than 10k umis we may be looking at a doublet 
#cells.to.keep=numi[paroi]<10000
#
##aggregate counts of cells with non-zero counts per gene 
#non.zero.cells=list()
#for(g in colnames(rMat)){
#    #g='YJR138W'
#    r=rMat[cells.to.keep,g]
#    a=aMat[cells.to.keep,g]
#    #efs=experiment.factor[paroi][numi[paroi]<5000]
#    ex=cbind(a,r)#[numi[paroi]<5000,] #a)
#    non.zero.cells[[g]]=sum(r+a>0)
#}
#
## some thresholding on how many informative cells we need to observe for each transcript
#genes.to.test=names(which(unlist(non.zero.cells)>64))
#
##beta binomial test for ASE
#bbin=mclapply(genes.to.test,
#    function(g,...){
#        r=rMat[cells.to.keep,g]
#        a=aMat[cells.to.keep,g]
#        #efs=experiment.factor[paroi][numi[paroi]<5000]
#        ex=cbind(a,r)#[numi[paroi]<5000,] #a)
#        #dcount[[g]]=colSums(ex)
#        print(g)
#        return(glmmTMB(ex~1,family=glmmTMB::betabinomial(link='logit'))) 
#    },
#    rMat,aMat,cells.to.keep,mc.cores=48)
#names(bbin)=genes.to.test
#bbin.llik=sapply(bbin, logLik)
#bbin.model.converged.genes=names(bbin.llik)[!is.na(bbin.llik)]
#
##why so slow??????
#lbb=lapply(bbin[bbin.model.converged.genes], function (x) coef(summary(x))$cond)
#rm(bbin)
#lbbdf=rbindlist(lapply(lbb, data.frame), idcol='gene')
#
#scASE_18=lbbdf
##write_csv(scASE_18, file='~/Dropbox/FileTransfer/scASE_18.csv')
#
#sig.genes=as.vector(unlist(scASE_18[,1][which(p.adjust(unlist(scASE_18[,5]), method='fdr')<.05)]))
#
#ncells=sum(cells.to.keep)
#cellID=factor(as.character(seq(1:ncells))) #length(r)
#setup.vars=list(
#        cellIDf=c(cellID,cellID),
#        geno=factor(c(rep('A',ncells),rep('B',ncells))),
#        of2=c(log(numi[paroi][cells.to.keep]),
#               log(numi[paroi][cells.to.keep]))
#        )
#
#nbin=mclapply(sig.genes, 
#   function(g, ... ){
#        print(g)
#        r=rMat[cells.to.keep,g]
#        a=aMat[cells.to.keep,g]
#        #efs=experiment.factor[paroi][numi[paroi]<5000]
#        y=c(r,a)#[numi
#        cellIDf=setup.vars$cellIDf
#        of2=setup.vars$of2
#        geno=setup.vars$geno
#        nbfit=glmmTMB(y~geno+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
#        if(!is.na(logLik(nbfit))){     nbin[[g]]=summary(nbfit)  } else {nbin[[g]]=NULL }
#   },
#  rMat,aMat,cells.to.keep,setup.vars,mc.cores=48)
#names(nbin)=sig.genes
#
#nbin.model.converged.genes=names(nbin)[sapply(nbin, function(x) !is.null(x))]
#scASE_18_nbin.mean=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$cond[2,])))
#scASE_18_nbin.varInt=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$disp[1,])))
#scASE_18_nbin.var=data.frame(gene=nbin.model.converged.genes,t(sapply(nbin[nbin.model.converged.genes], function(x) coef(x)$disp[2,])))
#
#scASE_18_nbin=left_join(scASE_18_nbin.mean, scASE_18_nbin.var, by="gene",suffix=c(".mean", ".var"))
#
#
#
#mElife_scASE_18=left_join(aseElife, scASE_18, by='gene')
#bbin_nbin=left_join(scASE_18, scASE_18_nbin,by='gene')
#
##compare bbin estimate to neg bin estimate 
#plot(unlist(bbin_nbin[,2]), unlist( bbin_nbin[,6]))
#abline(0,1)
#
##compare bbin fold change to elife paper ASE estimate
#plot(mElife_scASE_18$log2ASEFoldChange, log2(exp(mElife_scASE_18$Estimate)) )
#cor.test(mElife_scASE_18$log2ASEFoldChange[mElife_scASE_18$"Pr...z.."<.05],
#         log2(exp(mElife_scASE_18$Estimate[mElife_scASE_18$"Pr...z.."<.05] )) )
#cor.test(mElife_scASE_18$log2ASEFoldChange[mElife_scASE_18$"Pr...z.."<1],
#         log2(exp(mElife_scASE_18$Estimate[mElife_scASE_18$"Pr...z.."<1] )) )
#
#plot( unlist( bbin_nbin[,6]), unlist( bbin_nbin[,10]))
#
#hist(unlist( bbin_nbin[,13]), breaks=100)
#
#
#
#
#plot( unlist( bbin_nbin[,6]), unlist( bbin_nbin[,10]), col=(p.adjust(unlist( bbin_nbin[,13]), method='fdr')<.05)+1)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##tt
#
#
#
#
#
#
#
#
#
###capping at cells with less than 10k umis
##for(g in colnames(rMat)){
##    #g='YJR138W'
##    r=rMat[cells.to.keep,g]
##    a=aMat[cells.to.keep,g]
##    #efs=experiment.factor[paroi][numi[paroi]<5000]
##    ex=cbind(a,r)#[numi[paroi]<5000,] #a)
##    non.zero.cells[[g]]=sum(r+a>0)
##    #dcount[[g]]=colSums(ex)
##    print(g)
##    if(non.zero.cells[g]>64) {
##        bbin.fit=glmmTMB(ex~1,family=glmmTMB::betabinomial(link='logit'), control=glmmTMBControl(parallel=72))
##        print(bbin.fit)
##        if(!is.na(logLik(bbin.fit))){    bbin[[g]]=bbin.fit }
##    }
##}
#
##lbb=lapply(bbin, function (x) coef(summary(x))$cond)
##lbbdf=rbindlist(lapply(lbb, data.frame), idcol='gene')
#
#nbin=list()
##nbinCC=list()
#for(g in sig.genes ) { #colnames(rMat)){
#  #  if (g%in%colnames(rMat)){
#    print(g)
#    r=rMat[cells.to.keep,g]
#    a=aMat[cells.to.keep,g]
#    #efs=experiment.factor[paroi][numi[paroi]<5000]
#    y=c(r,a)#[numi
#    ncells=length(r)
#    cellID=factor(as.character(seq(1:ncells)))
#    cellIDf=c(cellID,cellID)
#    geno=factor(c(rep('A',ncells),rep('B',ncells)))
#    of2=c(log(numi[cells.to.keep]),
#            log(numi[cells.to.keep]))
#    if(non.zero.cells[g]>64) {
#        nbfit=glmmTMB(y~geno+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
#        if(!is.na(logLik(nbfit))){     nbin[[g]]=summary(nbfit)  }
#    #, control=glmmTMBControl(parallel=16))) #family=nbinom2))
#   }
#}
#
#scASE_18_nbin.mean=data.frame(gene=names(nbin),t(sapply(nbin, function(x) coef(x)$cond[2,])))
#scASE_18_nbin.varInt=data.frame(gene=names(nbin),t(sapply(nbin, function(x) coef(x)$disp[1,])))
#scASE_18_nbin.var=data.frame(gene=names(nbin),t(sapply(nbin, function(x) coef(x)$disp[2,])))
#
#
#scASE_18_nbin=left_join(scASE_18_nbin.mean, scASE_18_nbin.var, by="gene",suffix=c(".mean", ".var"))
#write_csv(scASE_18_nbin, file='~/Dropbox/FileTransfer/scASE_18_nbin.csv')
#
#
#cbind(scASE_18_nbin.mean, scASE_18_nbin.var)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#scASE_18=left_join(aseElife, lbbdf, by='gene')
#
#plot(scASE_18$log2ASEFoldChange[!is.na(scASE_18$z.value)], 
#     scASE_18$Estimate[!is.na(scASE_18$z.value)], ylim=c(-5,5))
#cor.test(scASE_18$log2ASEFoldChange[!is.na(scASE_18$z.value)], 
#         scASE_18$Estimate[!is.na(scASE_18$z.value)], method='spearman')
#
#dcount.df=do.call('rbind', dcount)
#plot(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#cor.test(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#
#y=c(r,a)
#geno=as.factor(c(rep('REF', length(r)), rep('ALT', length(a))))
#geno=relevel(geno, 'REF')
#cell=as.factor(c(seq_along(r),seq_along(r)))
#
#boxplot(c(r/totRcounts, a/totAcounts)~geno)
##of=log(numi[paroi])
##of2=c(of,of)
#of2=c(log(totRcounts), log(totAcounts))
#efs2=c(efs,efs)
#testNB=glmmTMB(y~geno+efs2+(1|cell)+offset(of2),family=nbinom2(link='log'))
#
##beta binomial
#ex=cbind(a,r) #a)
#testBB=glmmTMB(ex~efs,family=glmmTMB::betabinomial(link='logit'))
#lbbdf=rbindlist(lapply(lbb, data.frame), idcol='gene')
#
#aseElife=gdata::read.xls('/data/single_cell_eQTL/yeast/reference/elife-35471-data7-v2_ASE_results.xlsx')
#
#test=left_join(aseElife, lbbdf, by='gene')
#
#plot(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#cor.test(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#plot(test$log2ASEFoldChange, 
#     test$Estimate)
#plot(test$log2eQTLFoldChange[test$Pr...z..<.00001 ], 
#     test$Estimate[test$Pr...z..<.00001 ]) #, method='spearman')
#cor.test(test$log2eQTLFoldChange[test$Pr...z..<.0001 ], 
#     test$Estimate[test$Pr...z..<.0001 ], method='spearman')
#
#AC=getASEcounts(g.counts, sgd.genes, markerGR)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#p.assign.2=p.assign[to.keep]
#numiG.2=numiG[to.keep]
#numi.2=numi[to.keep]
#aC.2=aC[,to.keep]
#rC.2=rC[,to.keep]
#barcodes.2=barcodes[to.keep]
#counts.2=counts[,to.keep]
#
#
#p.assign.1=p.assign[to.keep]
#numiG.1=numiG[to.keep]
#numi.1=numi[to.keep]
#aC.1=aC[,to.keep]
#rC.1=rC[,to.keep]
#barcodes.1=barcodes[to.keep]
#counts.1=counts[,to.keep]
##library(session)
#
##save.session('/data/single_cell_eQTL/yeast/results/diploidASE_022520.RSession')
##restore.session('/data/single_cell_eQTL/yeast/results/diploidASE_022520.RSession')
#
#
##rownames(g.counts[[1]])
##m.granges=getMarkerGRanges(g.counts)
##m.granges$gcoord=gcoord.key[as.character(seqnames(m.granges))]+start(m.granges)
##m.granges$sname=colnames(vg)
#
##this above bit is to merge the two runs together 
#experiment.factor=factor(c(rep(1, length(p.assign.1)), rep(2,length(p.assign.2))))
#p.assign = c(p.assign.1, p.assign.2)
#numiG=c(numiG.1, numiG.2)
#numi=c(numi.1, numi.2)
#aC=cbind(aC.1, aC.2)
#rC=cbind(rC.1, rC.2)
#barcodes=c(barcodes.1, barcodes.2)
#counts=cbind(counts.1, counts.2)
#
###
#
#
#
## rejigger logic , 5/10/21
#
#cross.index=2
#crosses.to.parents[[cross.index]]
#phase_info=gts[match(paste0(variants[[1]], '_', variants[[2]]), rownames(gts)),]
#p1.ref= phase_info[, crosses.to.parents[[cross.index]][1]]=='0' 
#p1.ref[which(is.na(p1.ref))]=TRUE
#
#p1=rbind(rCs[p1.ref,],aCs[!p1.ref,])
#norder=c(which(p1.ref),which(!p1.ref))
#p1=p1[norder,]
#
##p1=p1[rownames(rC),]
#p2=rbind(rCs[!p1.ref,],aCs[p1.ref,])
#norder=c(which(!p1.ref),which(p1.ref))
#p2=p2[norder,]
#
##rownames(p1)=markerGR$name
#    #p2=p2[rownames(rC),]
#g.counts=list(p1,p2)
#
#m.granges=markerGR
#threads=16
#
# fO=findOverlaps(sgd.genes, m.granges)
#    mgsplit=split(m.granges$name[subjectHits(fO)], sgd.genes$Name[queryHits(fO)])
#    print('calc counts per gene, allele 1')
#    # could restructure these as sparse matrices if necessary
#    ref.ASEcounts=mcmapply( function(x) {
#            m=match(x, markerGR$name)                               
#            if(length(x)>1) {
#                colSums(g.counts[[1]][m,])
#            } else {
#                g.counts[[1]][m,]
#            }
#            #colSums(g.counts[[2]][x,]) 
#       },mgsplit,mc.cores=threads)
#
#    print('calc counts per gene, allele 2')
#    alt.ASEcounts=mcmapply( function(x) {
#            m=match(x, markerGR$name)                               
#            if(length(x)>1) {
#                colSums(g.counts[[2]][m,])
#            } else {
#                g.counts[[2]][m,]
#            }
#            #colSums(g.counts[[2]][x,]) 
#       }, mgsplit,mc.cores=threads)
# #       return(list(ref=ref.ASEcounts, alt=alt.ASEcounts))
##}
#
#
#                          
#                              
#
#
#
#
#
#
#
#
#
## extract each parents worth of data, flip counts, assign genes ...
#
#p1.ref= gts[,colnames(gts)[cross.index]]=='0' & (gts[,colnames(gts)[cross.index]]!=gts[,colnames(gts)[cross.index+1]]) 
#
#
#cross.index=2
#p1.ref=gts[,colnames(gts)[cross.index]]=='0'
#paroi=which(p.assign==cross.index)
#
#rCs=rC[,paroi]
#aCs=aC[,paroi]
#
#p1=rbind(rCs[p1.ref,],aCs[!p1.ref,])
#norder=c(which(p1.ref),which(!p1.ref))
#p1=p1[norder,]
#
##p1=p1[rownames(rC),]
#p2=rbind(rCs[!p1.ref,],aCs[p1.ref,])
#norder=c(which(!p1.ref),which(p1.ref))
#p2=p2[norder,]
#
##rownames(p1)=markerGR$name
#    #p2=p2[rownames(rC),]
#g.counts=list(p1,p2)
#
#threads=16
#m.granges=markerGR[gts[,colnames(gts)[cross.index]]=='0' & (gts[,colnames(gts)[cross.index]]!=gts[,colnames(gts)[cross.index+1]]) ,]
#
##getASEcounts=function(g.counts, sgd.genes, m.granges, threads=16) {
#    # can split this back out if we want to keep track of counts from individual variants     
#    fO=findOverlaps(sgd.genes, m.granges)
#    mgsplit=split(m.granges$name[subjectHits(fO)], sgd.genes$Name[queryHits(fO)])
#    print('calc counts per gene, allele 1')
#    # could restructure these as sparse matrices if necessary
#    ref.ASEcounts=mcmapply( function(x) {
#            m=match(x, markerGR$name)                               
#            if(length(x)>1) {
#                colSums(g.counts[[1]][m,])
#            } else {
#                g.counts[[1]][m,]
#            }
#            #colSums(g.counts[[2]][x,]) 
#       },mgsplit,mc.cores=threads)
#
#    print('calc counts per gene, allele 2')
#    alt.ASEcounts=mcmapply( function(x) {
#            m=match(x, markerGR$name)                               
#            if(length(x)>1) {
#                colSums(g.counts[[2]][m,])
#            } else {
#                g.counts[[2]][m,]
#            }
#            #colSums(g.counts[[2]][x,]) 
#       }, mgsplit,mc.cores=threads)
# #       return(list(ref=ref.ASEcounts, alt=alt.ASEcounts))
##}
#
#
#
##debug what's happening for genes with fewer than the expected number of counts 
#r.lengths=sapply(ref.ASEcounts, length)
#a.lengths=sapply(alt.ASEcounts, length)
#
#bgenes=unique(c( names(which(a.lengths< max(r.lengths))), names(which(r.lengths< max(r.lengths)))))
#ref.ASEcounts=ref.ASEcounts[!names(ref.ASEcounts)%in%bgenes]
#alt.ASEcounts=alt.ASEcounts[!names(alt.ASEcounts)%in%bgenes]
#rMat=do.call('cbind', ref.ASEcounts)
#aMat=do.call('cbind', alt.ASEcounts)
#
#plot(log2(rowSums(rMat)), log2(rowSums(aMat)))
#abline(0,1)
#
#plot(log2(colSums(rMat)), log2(colSums(aMat)))
#abline(0,1)
#
#
#totRcounts=rowSums(rMat)
#totAcounts=rowSums(aMat)
#
#totGcounts=colSums(rMat+aMat)
#
#plot(log2(rowSums(rMat)[numi[paroi]<5000]), log2(rowSums(aMat)[numi[paroi]<5000]))
#
#umiCap=5000
#
#bbin=list()
#dcount=list()
#non.zero.cells=list()
#for(g in colnames(rMat)){
#    #g='YJR138W'
#    r=rMat[numi[paroi]<5000,g]
#    a=aMat[numi[paroi]<5000,g]
#    y=c(a,r)#[numi
#    #efs=experiment.factor[paroi][numi[paroi]<5000]
#    ex=cbind(a,r)#[numi[paroi]<5000,] #a)
#
#    dcount[[g]]=colSums(ex)
#    print(g)
#    if(totGcounts[g]>64) {
#    non.zero.cells[[g]]=sum(r+a>0)
#
#    bbin[[g]]=glmmTMB(ex~1,family=glmmTMB::betabinomial(link='logit'))
#    }
#}
#lbb=lapply(bbin, function (x) coef(summary(x))$cond)
##dcount.df=do.call('rbind', dcount)
#lbbdf=rbindlist(lapply(lbb, data.frame), idcol='gene')
#lbbdf$non.zero.cells=unlist(non.zero.cells)
#test=left_join(aseElife, lbbdf, by='gene')
#
#
#waspASE=read.csv('/data/single_cell_eQTL/yeast/reference/allele_specific_merged_with_genes2.txt', sep='\t', header=F)
#
#
#rSumW=sapply(split(waspASE[,13], waspASE[,4]), sum)
#aSumW=sapply(split(waspASE[,14], waspASE[,4]), sum)
#
#wRatio=log2(aSumW/rSumW)
#wRatio=data.frame(gene=names(wRatio), ratio=as.vector(wRatio))
#test=left_join(wRatio, lbbdf, by='gene')
#
#test=test[is.finite(test$ratio),]
#sel=test$non.zero.cells>64 & test$Pr...z..<.001 #& test$sigDeciderBonferroni>0
#
#plot(test$ratio[sel],
#     test$Estimate[sel])
#cor.test(test$ratio[sel],
#     test$Estimate[sel] ) #,method='spearman')
#
#
#nbin=list()
#non.zero.cells=list()
#for(g in aseElife$gene[114:nrow(aseElife)]) { #colnames(rMat)){
#    if (g%in%colnames(rMat)){
#    print(g)
#    r=rMat[numi[paroi]<5000,g]
#    a=aMat[numi[paroi]<5000,g]
#    #efs=experiment.factor[paroi][numi[paroi]<5000]
#    y=c(a,r)#[numi
#    non.zero.cells[[g]]=sum(r+a>0)
#    ncells=length(r)
#    cellID=factor(as.character(seq(1:ncells)))
#
#    cellIDf=c(cellID,cellID)
#    geno=factor(c(rep('A',ncells),rep('B',ncells)))
#    of2=c(log(numi[numi[paroi]<5000]),
#            log(numi[numi[paroi]<5000]))
#
#    nbin[[g]]=summary(glmmTMB(y~geno+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))) #family=nbinom2))
#    }
#}
##restart at YBR041W
##YBR147W
#match('YBR162C', aseElife$gene)
#
#nbin.stat=t(sapply(nbin, function(x) coef(x)$cond[2,]))
#nbin.stat=data.frame(gene=rownames(nbin.stat), nbin.stat, non.zero.cells=unlist(non.zero.cells[rownames(nbin.stat)]))
#
#test=left_join(aseElife, nbin.stat, by='gene')
#plot(test$log2ASEFoldChange[test$non.zero.cells>50],  test$Estimate[test$non.zero.cells>50])
#cor.test(test$log2ASEFoldChange[test$non.zero.cells>100],  test$Estimate[test$non.zero.cells>100], method='spearman')
#
#
#
#test2=left_join(lbbdf, nbin.stat, by='gene')
#
#
#
#
#
#
#
#
#
#y=c(r,a)
#geno=as.factor(c(rep('REF', length(r)), rep('ALT', length(a))))
#geno=relevel(geno, 'REF')
#cell=as.factor(c(seq_along(r),seq_along(r)))
#
#boxplot(c(r/totRcounts, a/totAcounts)~geno)
##of=log(numi[paroi])
##of2=c(of,of)
#of2=c(log(totRcounts), log(totAcounts))
#efs2=c(efs,efs)
#testNB=glmmTMB(y~geno+efs2+(1|cell)+offset(of2),family=nbinom2(link='log'))
#
##beta binomial
#ex=cbind(a,r) #a)
#testBB=glmmTMB(ex~efs,family=glmmTMB::betabinomial(link='logit'))
#lbbdf=rbindlist(lapply(lbb, data.frame), idcol='gene')
#
#aseElife=gdata::read.xls('/data/single_cell_eQTL/yeast/reference/elife-35471-data7-v2_ASE_results.xlsx')
#
#test=left_join(aseElife, lbbdf, by='gene')
#
#plot(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#cor.test(test$log2ASEFoldChange[!is.na(test$z.value)], test$Estimate[!is.na(test$z.value)])
#
#plot(test$log2ASEFoldChange, 
#     test$Estimate)
#cor.test(test$log2eQTLFoldChange[test$Pr...z..<.001 ], 
#     test$Estimate[test$Pr...z..<.001 ], method='spearman')
#
#AC=getASEcounts(g.counts, sgd.genes, markerGR)
#
#
#library(MASS)
#library(VGAM)
#
#theta1=8
#theta2=.5
#ncells=500
#total.count.per.cell=5
#geno=factor(c(rep('A',ncells),rep('B',ncells)))
#cellID=factor(as.character(seq(1:ncells)))
#cellIDf=c(cellID,cellID)
#ycounts1=rnegbin(ncells, 1*total.count.per.cell,theta=theta1)
#ycounts2=rnegbin(ncells, 1.5*total.count.per.cell,theta=theta2)
#
#y=c(ycounts1,ycounts2)
#test=glmmTMB(y~geno+(1|cellIDf),family=nbinom2)
#test0=glmmTMB(y~geno+(1|cellIDf),family=nbinom2, dispformula=~geno)
#
#yy=cbind(ycounts2,ycounts1)
#test1=glmmTMB(yy~1, family=betabinomial(link='logit'))
## ASEdat has 3 columns: observed count for one allele, total coverage, "true" fold change
#
##https://www.biorxiv.org/content/10.1101/699074v2.full
#omDbetabinom=function(targetRho, ASEdat) {
#    -sum(dbetabinom(ASEdat$AlleleCount, ASEdat$Coverage, ASEdat$Prob, targetRho, log=TRUE))
#}
#ASEdat4OptimNoori <- data.frame(
#    NooriCountByGeneDownSampled[, 1],
#    rowSums(NooriCountByGeneDownSampled),
#    (jointASEData[,"NooriFC"] + jointASEData[,"AlbertFC"] / 2)
#)
#ASEdat4OptimNoori[,3] <- 2^ASEdat4OptimNoori[,3] / ( 1 + 2^ASEdat4OptimNoori[,3])
#names(ASEdat4OptimNoori) <- c("AlleleCount", "Coverage", "Prob")
#optimResultAlbert_all0.5 <- optim(0.5, omDbetabinom, ASEdat = ASEdat4OptimAll0.5_Albert, method='Brent', lower=1/1e8, upper=1)
#
#
#
#
#
#
#
#
#
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
#
#
#crosses.to.parents=list(
#     '375'=c("M22", "BYa"),
#     'A'  =c("BYa", "RMx"),
#     '376'=c("RMx", "YPS163a"),
#     'B'  =c("YPS163a", "YJM145x"),
#     '377'=c("YJM145x", "CLIB413a"),
#     '393'=c("CLIB413a", "YJM978x"),
#     '381'=c("YJM978x", "YJM454a"),
#    '3008'=c("YJM454a", "YPS1009x"),
#    '2999'=c("YPS1009x", "I14a"),
#    '3000'=c("I14a", "Y10x"),
#    '3001'=c("Y10x", "PW5a"),
#    '3049'=c("PW5a", "273614xa"),
#    '3003'=c("273614xa", "YJM981x"),
#    '3004'=c("YJM981x", "CBS2888a"),
#    '3043'=c("CBS2888a", "CLIB219x"),
#    '3028'=c("CLIB219x", "M22")
#    )
#
#
## for parallelizing the HMM
#ncores=8
#registerDoMC(cores=ncores)
## read in parental data ---------------------------------------------------------
#all.vcf=read.vcfR(paste0(reference.dir, 'parents.nostar.vcf'))
#gts=extract.gt(all.vcf)
#rm(all.vcf)
#
##par.vcf=read.vcfR(paste0(reference.dir, 'BYxRMxYPS163xYJM145.vcf'))
##gts=extract.gt(par.vcf)
#
##gts=gtsp[rownames(gts),]#c(2,9,11,16)]
#gts[is.na(gts)]="0"
##--------------------------------------------------------------------------------


