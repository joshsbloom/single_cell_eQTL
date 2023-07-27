
hotspotList=list()
for(set in names(sets)){
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))

   
    Gsub=segDataList$Gsub

    markerGR=getMarkerGRanges(list(t(Gsub)))


    cmin=sapply(sapply(split(markerGR, seqnames(markerGR)), start), min)
    cmin[!is.na(cmin)]=0
    cmax=sapply(sapply(split(markerGR, seqnames(markerGR)), start), max)

    cbin=data.frame(chr=names(cmin),
                    start=cmin,
                    end=cmax, 
                    strand="*")
    cbin=makeGRangesFromDataFrame(cbin)
    ###########################################################################################################3
    # analyze within each sub experiment 
    cbin50k=GenomicRanges::tile(cbin, width=50000)


   GP=readRDS(file=paste0(comb.out.dir,'/LOD_NB_combined_peaks.RDS'))
  
    # plotEQTL(GP[GP$LOD>4,], titleme='combined', CI=F)
   # cPf=GP[GP$LOD>4 & GP$chr!=GP$tchr,]
   cPf=GP[GP$FDR<.05 & GP$chr!=GP$tchr,]
    
   bin.table=makeBinTable(cPf, cbin50k)
   cH=plotHotspot2(bin.table,titleme=set)

   cPf=getBinAssignment(cPf, cbin50k)
   sig.hp=qpois(1-(.05/length(bin.table$pos)),ceiling(mean(bin.table$count)))+1

   sig.hp.names=table(cPf$bin)>sig.hp
   cPf$in.hotspot=sig.hp.names[cPf$bin]
   attr(cPf, 'threshold')=sig.hp

    plist=list()
    plist[['combined']]=cPf

    hlist=list()
    for(cc in cycle.cats) {
         GP=readRDS(paste0(comb.out.dir,'/LOD_NB_cell_cycle_peaks_', cc,'.RDS'))
         #    cPf=GP[GP$LOD>4 & GP$chr!=GP$tchr,]
         cPf=GP[GP$FDR<.05 & GP$chr!=GP$tchr,]

         bin.table=makeBinTable(cPf, cbin50k)
         cPf=getBinAssignment(cPf, cbin50k)
         sig.hp=qpois(1-(.05/length(bin.table$pos)),ceiling(mean(bin.table$count)))+1
         sig.hp.names=table(cPf$bin)>sig.hp
         cPf$in.hotspot=sig.hp.names[cPf$bin]
            attr(cPf, 'threshold')=sig.hp

         hlist[[cc]]=plotHotspot2(bin.table, titleme=cc)
         plist[[cc]]=cPf
     }
    
    ggpubr::ggarrange(cH, hlist[[1]], hlist[[2]], hlist[[3]], hlist[[4]], hlist[[5]],ncol=1)
        ggsave(file=paste0(comb.out.dir,'/LOD_NB_hotspots.png'), width=8, height=13)
    saveRDS(plist, file=paste0(comb.out.dir,'/hotspot_peaks.RDS'))
    hotspotList[[set]]=plist
}

for(set in names(sets)){
    print(set)
   comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    tTL=data.frame(readRDS(file=paste0(comb.out.dir,'/hotspot_peaks.RDS'))$combined)
    ucontrasts=combn(cycle.cats,2)
    ucontrastsB=matrix(gsub(':', '.', paste0('Beta_', ucontrasts)), nrow=2)
    ucontrastsSE=matrix(gsub(':', '.', paste0('SE_', ucontrasts)), nrow=2)

    tTL.contrasts=list()
    for (i in 1:ncol(ucontrasts)) {
        coi=paste(ucontrasts[1,i], ucontrasts[2,i])

        B1=tTL[,match(ucontrastsB[1,i], names(tTL))]
        B2=tTL[,match(ucontrastsB[2,i], names(tTL))]

        SE1=tTL[,match(ucontrastsSE[1,i], names(tTL))]
        SE2=tTL[,match(ucontrastsSE[2,i], names(tTL))]

        Zstat=(B1-B2)/sqrt(SE1^2+SE2^2)
        pval=2*pnorm(-abs(Zstat),lower.tail=T)
        tTL.contrasts[[coi]]=data.frame(chrom=tTL$chrom, transcript=tTL$transcript, peak.marker=tTL$peak.marker,
                                        Zstat=Zstat, pval=pval)
    }
    sum(p.adjust(unlist(lapply(tTL.contrasts, function(x) x$pval)), method='fdr')<.05)
    sort(unique(unlist(lapply(tTL.contrasts, function(x) x$transcript))[which(p.adjust(unlist(lapply(tTL.contrasts, function(x) x$pval)), method='fdr')<.05)]) )
    
    saveRDS(tTL.contrasts, file=paste0(comb.out.dir,'/trans_CC_test_CombinedResultsContrasts.RDS'))
}

visualize_CC_trans=function(tTL, tTL.contrasts) {
    top.sig=sort(unique(unlist(lapply(tTL.contrasts, function(x) paste0(x$transcript, '__', x$peak.marker)))[which(p.adjust(unlist(lapply(tTL.contrasts, function(x) x$pval)), method='fdr')<.001)]) )

    tTL.B=pivot_longer(tTL, col=starts_with("Beta_"), names_to="cc",  values_to="Beta")
    tTL.SE=pivot_longer(tTL, col=starts_with("SE_"), names_to="cc",  values_to="SE")
    tTL.l=tTL.B
    tTL.l$SE=tTL.SE$SE

    #tTL.l=tTL.l[,-match('cc2',names(tTL.l))]
    tTL.l$cc=gsub('Beta_' ,'',tTL.l$cc)
    tTL.l$cc=gsub('\\.' ,'/',tTL.l$cc)
    tTL.l$cc=factor(tTL.l$cc, levels=c('M/G1', 'G1', 'G1/S', 'S', 'G2/M'))
  
    tTL.l$tp=paste0(tTL.l$transcript, '__', tTL.l$peak.marker)
    tTL.l %>% filter(tp %in% top.sig) %>% ggplot() + geom_bar(aes(x=cc, y=Beta), stat='identity')+
            geom_errorbar(aes(x=cc, ymin=Beta-(1.96*SE), ymax=Beta+(1.96*SE)), width=.4)+facet_wrap(~transcript+peak.marker)+
             ggtitle(paste(set, '_____CC x Geno trans-eQTL Betas'))
}


#        ggsave(file=(paste0(results.base.dir,experiment.name,'/',dip,'-bbin_sigASE_CC.png')))


#lapply(hotspotList, function(x) { unique(x$combined$bin[x$combined$in.hotspot]) })

#test for cell_cycle interactions within hotspot


