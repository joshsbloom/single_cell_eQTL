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
library(mgcv)
library(Rfast2)
library(fastglm)
library(MASS)
library(irlba)
library(nebula)
library(dplyr)

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

source(paste0(code.dir, 'ASE_fxs.R'))

plotEQTL=function(cPf, titleme='',CI=T) {
  a=ggplot(cPf,aes(x=pos,y=tpos)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
        geom_point(size=.5)+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
        facet_grid(tchr~chr,scales="free", space='free')+
        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))+
        ggtitle(titleme)
    if(CI) {
        a=a+geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) 
    }
    return(a)
}




chroms=paste0('chr', as.roman(1:16)) 

sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))
sgd.genes=getSGD_GeneIntervals(reference.dir)
sgd.genes$chromosome_name=seqnames(sgd.genes)
sgd.genes$start_position=start(sgd.genes)
bulk_peaks=data.table::rbindlist(peaks.per.gene.augmented)

bulk_peaks$chr=tstrsplit(bulk_peaks$pmarker, ':')[[1]]
bulk_peaks$pos=as.numeric(tstrsplit(bulk_peaks$pmarker, ':|_')[[2]])
mbg=(match(bulk_peaks$gene,sgd.genes$gene_id))
bulk_peaks=bulk_peaks[!is.na(mbg),]
bulk_peaks$tchr=as.character(sgd.genes$chromosome_name[as.vector(match(bulk_peaks$gene,sgd.genes$gene_id))])
bulk_peaks$tpos=as.numeric(sgd.genes$start_position[as.vector(match(bulk_peaks$gene,sgd.genes$gene_id))])
    bulk_peaks$tchr=factor(bulk_peaks$tchr,levels=rev(chroms)) #c("X","V","IV","III","II","I"))
    bulk_peaks$chr=factor(bulk_peaks$chr,levels=chroms)    


comb.out.dir="/data/single_cell_eQTL/yeast/results/combined/A/"
sc_peaks=readRDS(file=paste0(comb.out.dir,'/LOD_NB_combined_peaks.RDS'))



a=plotEQTL(sc_peaks[sc_peaks$FDR<.05,], CI=F)

b=plotEQTL(bulk_peaks[bulk_peaks$var.exp>.04,], CI=F)

ggpubr::ggarrange(a,b, nrow=2, ncol=1)
  ggsave( '~/A_map_comp.png', width=20, height=10)



sc_peaks
        chrom transcript             peak.marker   CI.l   CI.r     LOD      FDR
    1:   chrI    YAL055W      chrI_41802_G_A_452   1483  42527   4.579 0.537805
    2:   chrI    YAL054C      chrI_48115_G_C_499   1483 138257   2.163 0.998981
    3:   chrI    YAL053W      chrI_48115_G_C_499   1483 228756   1.334 1.000000
    4:   chrI    YAL051W    chrI_228756_A_G_1542   1483 228756   2.271 0.998932
    5:   chrI    YAL049C      chrI_48115_G_C_499  48115  48115 156.841 0.000000
   ---                                                                         
65916: chrXVI    YPR193C chrXVI_943160_G_T_49701 920598 943160  10.314 0.001989
65917: chrXVI    YPR196W chrXVI_943160_G_T_49701 943160 943160  10.155 0.001989
65918: chrXVI    YPR198W chrXVI_718898_T_C_49323 694515 779281   3.740 0.897820
65919: chrXVI    YPR199C chrXVI_228950_G_A_48017 199138 916830   3.093 0.993631
65920: chrXVI    YPR200C chrXVI_943160_G_T_49701 173535 943160   2.787 0.997766
          chr    pos   tchr   tpos
    1:   chrI  41802   chrI  42177
    2:   chrI  48115   chrI  42761
    3:   chrI  48115   chrI  45899
    4:   chrI 228756   chrI  48564
    5:   chrI  48115   chrI  51737
   ---                            
65916: chrXVI 943160 chrXVI 922789
65917: chrXVI 943160 chrXVI 931376
65918: chrXVI 718898 chrXVI 934034
65919: chrXVI 228950 chrXVI 937982
65920: chrXVI 943160 chrXVI 939159

