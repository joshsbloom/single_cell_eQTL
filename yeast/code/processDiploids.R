#library(Matrix)
#library(data.table)
#library(qtl)
#library(qtl2)
#library(parallel)
#library(Rfast)
#library(seqinr)
#library(GenomicRanges)
##library("BSgenome.Scerevisiae.UCSC.sacCer3")
#library(foreach)
#library(doMC)
#library(ggplot2)
#library(vcfR)
#library(monocle)
#library(matrixStats)
#library(spam)
#
library(glmmTMB)
library(broom.mixed)
library(foreach)
library(doMC)
library(dplyr)

code.dir='/data/single_cell_eQTL/yeast/code/'
reference.dir='/data/single_cell_eQTL/yeast/reference/'
data.base.dir='/data/single_cell_eQTL/yeast/processed/'
results.base.dir='/data/single_cell_eQTL/yeast/results/'

source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))

#load custom functions for processing ASE data 
source(paste0(code.dir, 'ASE_fxs.R'))

sgd.genes=getSGD_GeneIntervals(reference.dir)

experiments=list(
    '12_Group_1_diploids_3004_2444_3051_5k_Feb_21'=names(crosses.to.parents)[c(2,4,14)],
    '13_Group1_diploids_3004_2444_3051_7point5K_Feb_21'=names(crosses.to.parents)[c(2,4,14)],
    '14_GP1_3004_3051_274_375_May10'=names(crosses.to.parents)[c(1,2,4,14)],
    '15_GP2_376_377_393_3008_May10'=names(crosses.to.parents)[c(3,5,6,8)],
    '16_GP3_2999_3000_3001_3049_May10'=names(crosses.to.parents)[c(9,10,11,12)],
    '17_GP4_3003_3043_3028_381_May10'=names(crosses.to.parents)[c(7,13,15,16)],
    '18_3051_May10'=names(crosses.to.parents)[2])

#good crosses A, B, 3004, 3008, 2999, 3001, 3003, 3028

good.dips=list(
    '12_Group_1_diploids_3004_2444_3051_5k_Feb_21'=names(crosses.to.parents)[c(2,4,14)],
    '13_Group1_diploids_3004_2444_3051_7point5K_Feb_21'=names(crosses.to.parents)[c(2,4,14)],
    '14_GP1_3004_3051_274_375_May10'=names(crosses.to.parents)[c(2,4,14)],
    '15_GP2_376_377_393_3008_May10'=names(crosses.to.parents)[c(8)],
    '16_GP3_2999_3000_3001_3049_May10'=names(crosses.to.parents)[c(9,11)],
    '17_GP4_3003_3043_3028_381_May10'=names(crosses.to.parents)[c(13,16)],
    '18_3051_May10'=names(crosses.to.parents)[2])


for(experiment.name in names(experiments)){
    input.diploids=experiments[[experiment.name]]
    print(input.diploids)
    #for output 
    dir.create(paste0(results.base.dir,experiment.name))
    data.dir=paste0(data.base.dir, experiment.name, '/')

    #vname=paste0(variants[[1]],'_', variants[[2]])
    ase.Data=buildASE_data(data.dir, experiment.name)
    saveRDS(ase.Data, paste0(results.base.dir,experiment.name,'/','aseData.RDS'))

    dip.Assignments=getDipAssignments(ase.Data, input.diploids, crosses.to.parents, ncores=8)
    dip.specificCounts=countDiploidSpecificVariants(ase.Data, input.diploids, crosses.to.parents)
    #overwrite dip.Assignments to include likelihood, counts, and umap
    dip.Assignments=dplyr::left_join(dplyr::left_join(ase.Data$um, dip.Assignments, by='barcode'), dip.specificCounts)
    saveRDS(dip.Assignments, paste0(results.base.dir,experiment.name,'/','diploid_assignments.RDS'))

    for(dip in input.diploids) {
        #dip='A'
        phasedCounts=getPhasedCountsPerTranscript(ase.Data, dip.Assignments, dip, sgd.genes, crosses.to.parents) 
        saveRDS(phasedCounts, file=paste0(results.base.dir,experiment.name,'/','phasedCounts_',dip,'.RDS'))

        bbin.model.results=doBetaBinomialTest(dip, phasedCounts,dip.Assignments, ase.Data)
        saveRDS(bbin.model.results, file=paste0(results.base.dir,experiment.name,'/','bbin_',dip,'.RDS'))
        
        nbin.model.results=doNbinTest(dip,phasedCounts,dip.Assignments,ase.Data)
        saveRDS(nbin.model.results, file=paste0(results.base.dir,experiment.name,'/','nbin_',dip,'.RDS'))
    }
}

aseElife=gdata::read.xls('/data/single_cell_eQTL/yeast/reference/elife-35471-data7-v2_ASE_results.xlsx')

experiment.name=names(experiments)[4]

library(ggplot2)
library(ggrepel)
library(viridis)
nUMI.thresh=10000

for(experiment.name in names(experiments)){
    input.diploids=experiments[[experiment.name]]
    dip.Assignments=readRDS(file=paste0(results.base.dir,experiment.name,'/','diploid_assignments.RDS'))
    ase.Data=readRDS(paste0(results.base.dir,experiment.name,'/','aseData.RDS'))
        
    print(input.diploids)
    for(dip in input.diploids) {
        phasedCounts=readRDS(file=paste0(results.base.dir,experiment.name,'/','phasedCounts_',dip,'.RDS'))
     #   plot(log2(rowSums(phasedCounts[[1]])), log2(rowSums(phasedCounts[[2]])))
             classified.cells=dip.Assignments$diploid_name==dip
        selected.barcodes=dip.Assignments$barcode[dip.Assignments$diploid_name==dip]
        cells.to.keep=ase.Data$numi[classified.cells]<nUMI.thresh 
        p1=phasedCounts[[1]][cells.to.keep,]
        p2=phasedCounts[[2]][cells.to.keep,]
        gInfoTotals=p1+p2 #phasedCounts[[1]][cells.to.keep,]+phasedCounts[[2]][cells.to.keep,]
        informativeCellsPerTranscript=colSums(gInfoTotals>0)
        paste0(results.base.dir,experiment.name,'/','diploid_assignments.RDS')
        infoCells=data.frame(gene=names(informativeCellsPerTranscript), genoInfoCells=as.vector(informativeCellsPerTranscript), stringsAsFactors=F)
        

        parent1=crosses.to.parents[[dip]][1]
        parent2=crosses.to.parents[[dip]][2]


        png(paste0(results.base.dir,experiment.name,'/',dip, '-countsPerTranscript.png'), width=768,height=768)
        plot(log2(colSums(phasedCounts[[1]])), log2(colSums(phasedCounts[[2]])),
             main=paste(experiment.name, dip),
             sub='geno informative counts per transcript',
             xlab=paste('log2', parent1),
             ylab=paste('log2', parent2))
        abline(0,1)
        dev.off()

        gCounts=sapply(phasedCounts, colSums)
        png(paste0(results.base.dir,experiment.name,'/',dip,'-countsPerTranscript_GenomePos.png'), width=1024,height=512)

        plot(sgd.genes$gcoord[match(rownames(gCounts),sgd.genes$gene_id)], 
             log2(gCounts[,1]/gCounts[,2]),
             xlab='genomic coordinate per each transcript',
             ylab=paste0('log2(', parent1, '/', parent2, ')'),
             main=paste(experiment.name, dip)   )
        abline(v=gcoord.key, lty=2, col='lightblue')
        abline(h=0, col='black')
        dev.off()
    }
}




for(experiment.name in names(experiments)[1:6]){
    input.diploids=good.dips[[experiment.name]]
    dip.Assignments=readRDS(file=paste0(results.base.dir,experiment.name,'/','diploid_assignments.RDS'))
    ase.Data=readRDS(paste0(results.base.dir,experiment.name,'/','aseData.RDS'))
        
    print(input.diploids)
    for(dip in input.diploids) {
        phasedCounts=readRDS(file=paste0(results.base.dir,experiment.name,'/','phasedCounts_',dip,'.RDS'))
     #   plot(log2(rowSums(phasedCounts[[1]])), log2(rowSums(phasedCounts[[2]])))
        classified.cells=dip.Assignments$diploid_name==dip
        selected.barcodes=dip.Assignments$barcode[dip.Assignments$diploid_name==dip]
        cells.to.keep=ase.Data$numi[classified.cells]<nUMI.thresh 
        p1=phasedCounts[[1]][cells.to.keep,]
        p2=phasedCounts[[2]][cells.to.keep,]
        gInfoTotals=p1+p2 #phasedCounts[[1]][cells.to.keep,]+phasedCounts[[2]][cells.to.keep,]
        informativeCellsPerTranscript=colSums(gInfoTotals>0)
        paste0(results.base.dir,experiment.name,'/','diploid_assignments.RDS')
        infoCells=data.frame(gene=names(informativeCellsPerTranscript), genoInfoCells=as.vector(informativeCellsPerTranscript), stringsAsFactors=F)
        

        parent1=crosses.to.parents[[dip]][1]
        parent2=crosses.to.parents[[dip]][2]

   
       bbin.model.results=readRDS(paste0(results.base.dir,experiment.name,'/','bbin_',dip,'.RDS'))
       nbin.model.results=readRDS(paste0(results.base.dir,experiment.name,'/','nbin_',dip,'.RDS'))
       #left_join(aseElife, nbin.model.results %>% filter(term=='genoB' & component=='cond'), by='gene')

       nbin.model.resultsM=left_join(nbin.model.results %>% filter(term=='genoB' & component=='cond'),
                                      nbin.model.results %>% filter(term=='genoB' & component=='disp'), by='gene', suffix=c('.cond', '.disp'))

       # if likelihood ratio test
        nbin.model.resultsM=nbin.model.resultsM %>% mutate(p.llrt.disp=pchisq(-2*(logLik.red.disp-logLik.disp),1,lower.tail=F))
       # ggplot(nbin.model.resultsM, aes(x=-log10(p.value.disp), y=-log10(p.llrt.disp)))+geom_point()+theme_bw()
        

        all_models=left_join(bbin.model.results %>% filter(term=='(Intercept)'), nbin.model.resultsM, by='gene')
        all_models=left_join(all_models, infoCells, by='gene')

        ggplot(all_models, aes(x=estimate, y=-log10(p.value), col=log2(genoInfoCells)))+geom_point()+xlim(c=-10,10)+
             scale_colour_viridis() +
               geom_label_repel(aes(x=estimate,
                                    y=-log10(p.value),
                               label=ifelse(p.adjust(p.value, method='fdr')<.05,gene,'')))+theme_bw()+
                               xlab(paste('Mean effect beta-binomial', parent2, 'vs', parent1))+
             ggtitle(paste(experiment.name, dip))
        ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-bbin_sigASE.png')))



        #replace with logLik.disp
        a=ggplot(all_models %>% filter(!is.na(logLik)) %>% mutate(sig.disp=p.adjust(p.value.disp,'fdr')<.05) ,
               aes(x=estimate.disp, y=estimate.cond,  #paste(p.adjust(p.value.disp,'fdr')<.05),
               col=log2(genoInfoCells))) + #log2(genoInfoCells))) +
            geom_point()+
          scale_colour_viridis() + # name='lodispersion q<.05') +
          #scale_colour_manual(name='wald dispersion q<.05', values=c('black','red')) +
            xlab(paste('Dispersion effect negative-binomial', parent2, 'vs', parent1)) + 
            ylab(paste('Mean effect negative-binomial', parent2, 'vs', parent1)) +  theme_bw()+
            geom_label_repel(aes(x=estimate.disp,
                               y=estimate.cond,
                               label=ifelse(sig.disp,gene,'')))+theme_bw()+
             ggtitle(paste(experiment.name, dip))
         b=ggplot(all_models %>% filter(!is.na(logLik)), aes(x=p.value.disp))+geom_histogram()+ xlab('(p) wald disp effect')+theme_bw()
         ggarrange(a,b, nrow=2, heights=c(1,.25))
         ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-nbin_sigDisp_ASE.png')))
    }
}



#lower in BY
to.test=nbin.model.resultsM[which(-log10(nbin.model.resultsM$p.llrt.disp)>50),]
g='YLR350W'
g='YAR019C'
g='YBR031W'
g='YNL336W'
g='YBR012C'
g='YJL201W'
g='YKL160W'
g='YBR043C'

#YDL207W
for(g in to.test$gene) {
    print(g)
    aseElife[aseElife$gene==g,]
    all_models[all_models$gene==g,]
    plot( jitter(log10(p1[,g]+1)), jitter(log10(p2[,g]+1))
         , xlab=parent1, ylab=parent2)
    readline()

   plot( (log10((p1[,g]+1)/rowSums(p1+p2))), (log10((p2[,g]+1)/rowSums(p1+p2)))
    , xlab=parent1, ylab=parent2)
    readline()

}


    y=c(log10(p1[,g]+1/rowSums(p1+p2)),log10(p2[,g]+1/rowSums(p1+p2)))
    x=rep(c(parent1,parent2),each=nrow(p1))
    #readline()

    ggb=data.frame(geno=x, ncounts=y)
    a=ggplot(ggb,aes(geno,ncounts))+geom_quasirandom()
      plot(a)
}
    
       bbin.model.results=readRDS(paste0(results.base.dir,experiment.name,'/','bbin_',dip,'.RDS'))
       nbin.model.results=readRDS(paste0(results.base.dir,experiment.name,'/','nbin_',dip,'.RDS'))
       #left_join(aseElife, nbin.model.results %>% filter(term=='genoB' & component=='cond'), by='gene')

        nbin.model.resultsM=left_join(nbin.model.results %>% filter(term=='genoB' & component=='cond'),
                                      nbin.model.results %>% filter(term=='genoB' & component=='disp'), by='gene', suffix=c('.cond', '.disp'))

        # if likelihood ratio test
        nbin.model.resultsM=nbin.model.resultsM %>% mutate(p.llrt.disp=pchisq(-2*(logLik.red.disp-logLik.disp),1,lower.tail=F))

        a=ggplot(nbin.model.resultsM, aes(x=-log10(p.value.disp), y=-log10(p.llrt.disp)))+geom_point()+theme_bw()+
            xlab('-log10(p) wald') +
            ylab('-log10(p) lrt') +ggtitle('Dispersion test comparison')

        b=ggplot(nbin.model.resultsM, aes(x=p.value.disp))+geom_histogram()+ xlab('(p) wald')
        c=ggplot(nbin.model.resultsM, aes(x=p.llrt.disp))+geom_histogram()+xlab('(p) lrt')
          grid.arrange(a,b,c, ncol=2, nrow=2, layout_matrix=rbind(c(1,2),c(1,3)))+ggtitle(paste(experiment.name, dip))
        ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-compWaldLRT.png')))
        

        all_models=left_join(bbin.model.results %>% filter(term=='(Intercept)'), nbin.model.resultsM, by='gene')
        all_models=left_join(all_models, infoCells, by='gene')

        ggplot(all_models %>% filter(!is.na(logLik.disp)) %>% mutate(sig.disp=p.adjust(p.value.disp,'fdr')<.05) %>%
               aes(x=estimate.disp, y=estimate.cond,  #paste(p.adjust(p.value.disp,'fdr')<.05),
               col=sig.disp)) + #log2(genoInfoCells))) +
            geom_point()+
         #scale_colour_viridis() + # name='lodispersion q<.05') +
          scale_colour_manual(name='wald dispersion q<.05', values=c('black','red')) +
            xlab(paste('Dispersion effect negative-binomial', parent2, 'vs', parent1)) + 
            ylab(paste('Mean effect negative-binomial', parent2, 'vs', parent1)) +  theme_bw()+
            geom_label_repel(aes(x=estimate.disp,
                               y=estimate.cond,
                               label=ifelse(sig.disp,gene,'')))+theme_bw()+
             ggtitle(paste(experiment.name, dip))
      
         ggplot(all_models %>% filter(!is.na(logLik.disp)), 
                aes(x=estimate.cond, y=-log10(p.value.cond)))+geom_point()
       

           doplot=F
        if(experiment.name == '18_3051_May10' & doplot==T) {
        #plot(all_models$estimate[!is.na(all_models$logLik.disp)], 
        #     all_models$estimate.cond[!is.na(all_models$logLik.disp)], xlim=c(-5,5), ylim=c(-5,5),
        #     xlab=paste('bbin log(', parent2, '/', parent1, ')'),
        #     ylab=paste('nbin log(', parent2, '/', parent1, ')'))
        #abline(0,1)
     
         #compare negbin vs bbin
        ggplot(all_models %>% filter(!is.na(logLik.disp)),
               aes(x=estimate, y=estimate.red.cond, col=paste(p.adjust(p.value,'fdr')<.05)))+
            geom_point()+
          scale_colour_manual(name='beta binomial q<.05', values=c('black','red')) +
            xlab(paste('Mean effect beta-binomial', parent2, 'vs', parent1)) + 
            ylab(paste('Mean effect negative-binomial', parent2, 'vs', parent1)) + 
            #geom_label_repel(aes(x=estimate,
            #                   y=estimate.cond,
            #                   label=ifelse(p.adjust(p.value, 'fdr')<.05,gene,'')))+
           geom_abline(slope=1, intercept=0)+
            theme_bw()
        ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-bbin_v_negbin.png')))


        #plot(all_models$estimate.disp[!is.na(all_models$logLik.disp)], 
        #     all_models$estimate.cond[!is.na(all_models$logLik.disp)], 
        #     xlab=paste('Dispersion effect negative-binomial', parent2, 'vs', parent1),
        #     ylab=paste('Mean effect negative-binomial', parent2, 'vs', parent1),
        #     col=(p.adjust(all_models$p.value.disp[!is.na(all_models$logLik.disp)],'fdr')<.05)+1 )
        #abline(0,1)

        a=ggplot(all_models %>% filter(!is.na(logLik.disp)) %>% mutate(sig.disp=p.adjust(p.value.disp,'fdr')<.05) %>%
                 mutate(sig.displl=p.adjust(p.llrt.disp, 'fdr')<.05),
               aes(x=estimate.disp, y=estimate.cond,  #paste(p.adjust(p.value.disp,'fdr')<.05),
               col=sig.disp)) + #log2(genoInfoCells))) +
            geom_point()+
        #scale_colour_viridis() + # name='lodispersion q<.05') +
          scale_colour_manual(name='wald dispersion q<.05', values=c('black','red')) +
            xlab(paste('Dispersion effect negative-binomial', parent2, 'vs', parent1)) + 
            ylab(paste('Mean effect negative-binomial', parent2, 'vs', parent1)) +  theme_bw()+
          #  geom_label_repel(aes(x=estimate.disp,
          #                     y=estimate.cond,
          #                     label=ifelse(sig.disp,gene,'')))+theme_bw()+
             ggtitle(paste(experiment.name, dip))

        b=ggplot(all_models %>% filter(!is.na(logLik.disp)) %>% mutate(sig.disp=p.adjust(p.value.disp,'fdr')<.05) %>%
                 mutate(sig.displl=p.adjust(ifelse(p.llrt.disp==0, 2e-16,p.llrt.disp), 'fdr')<.05),
               aes(x=estimate.disp, y=estimate.cond,  #paste(p.adjust(p.value.disp,'fdr')<.05),
               col=sig.displl)) + #log2(genoInfoCells))) +
            geom_point()+
        #scale_colour_viridis() + # name='lodispersion q<.05') +
          scale_colour_manual(name='llrt dispersion q<.05', values=c('black','red')) +
            xlab(paste('Dispersion effect negative-binomial', parent2, 'vs', parent1)) + 
            ylab(paste('Mean effect negative-binomial', parent2, 'vs', parent1)) +  theme_bw()+
          #  geom_label_repel(aes(x=estimate.disp,
          #                     y=estimate.cond,
          #                     label=ifelse(sig.displl,gene,'')))+theme_bw()+
             ggtitle(paste(experiment.name, dip))



        c= ggplot(all_models %>% filter(!is.na(logLik.disp)),
               aes(x=estimate.disp, y=estimate.cond, col=log2(genoInfoCells)))+
            geom_point()+
        scale_colour_viridis() + # name='lodispersion q<.05') +
            xlab('Dispersion effect RM v BY') + 
            ylab('Mean effect RM v BY') + theme_bw()+
            #geom_label_repel(aes(x=estimate.disp,
            #                   y=estimate.cond,
            #                   label=ifelse(p.adjust(p.value.disp,'fdr')<.05,gene,'')))+theme_bw()+
             ggtitle(paste(experiment.name, dip))

    ggarrange(a,b,c, nrow=1)
    ggsave(paste0(results.base.dir,experiment.name,'/',dip,'-mean_v_disp_comp.png'))



        aseElifeMerged=left_join(aseElife, bbin.model.results %>% filter(term=='(Intercept)'), by='gene')
        aseElifeMerged=aseElifeMerged[!is.na(aseElifeMerged$p.value),]
        ggplot(aseElifeMerged,
               aes(x=log2(exp(estimate)), y=log2ASEFoldChange, col=paste(p.adjust(p.value,'fdr')<.05, sigDeciderFDR>0)))+
            geom_point()+
        scale_colour_discrete(name='| sig sc | sig bulk |') +
            xlab('log2(bbin estimate sc)') + 
            ylab('log2ASEFoldChange bulk') + 
            geom_label_repel(aes(x=log2(exp(estimate)),
                               y=log2ASEFoldChange,
                               label=ifelse(p.adjust(p.value,'fdr')<.05,gene,'')))+theme_bw()+
        ggtitle(paste(experiment.name, dip, 'r = ', 
                      round(cor(log2(exp(aseElifeMerged$estimate)),aseElifeMerged$log2ASEFoldChange),2)))
        ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-compBulk.png')))

        }
    }
}  




#bbin.model.results %>% filter(term=='(Intercept)')

bbin.model.results$p.value[bbin.model.results$p.value==0]=2e-16

plot(log2(exp(test$estimate)),test$log2ASEFoldChange)

test=left_join(aseElife, nbin.model.results %>% filter(term=='genoB' & component=='cond'), by='gene')
plot(log2(exp(test$estimate)),test$log2ASEFoldChange)

test=left_join(aseElife, nbin.model.results %>% filter(term=='genoB' & component=='disp'), by='gene')
plot(log2(exp(test$estimate)),test$log2ASEFoldChange)

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


