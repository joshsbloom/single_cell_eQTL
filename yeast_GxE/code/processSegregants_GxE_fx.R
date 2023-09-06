makeExpResults=function(base.dir,sets,experiments,data.dirs, chroms, state.annot.file, cc.bad.vars=c('HAPLOID','HAPLOIDS', 'MATALPHA', '?')){
 
    exp.results=list()

    for( ee in  sets[[set]] ) { 
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
        classification=readRDS(paste0(results.dir, 'cellFilter.RDS'))
        het.cells=!classification

          tmp=as_matrix(g.counts$ref.counts[,-het.cells])
          rsum=Rfast::rowsums(tmp)
          rm(tmp)
          tmp=as_matrix(g.counts$alt.counts[,-het.cells])
          asum=Rfast::rowsums(tmp)
          rm(tmp)
          #rsum=apply(g.counts$ref.counts[,-het.cells], 1, sum)
            #asum=apply(g.counts$alt.counts[,-het.cells],1, sum)
         af=(rsum/(rsum+asum))
         names(af)=rownames(g.counts$ref)

        #rsum=apply(g.counts$ref.counts,1, sum)
        #asum=apply(g.counts$alt.counts,1, sum)
        #af=(rsum/(rsum+asum))
        #af[(rsum+asum)<10]=NA
        af=data.frame(chr=tstrsplit(names(af), '_')[[1]], pos=as.numeric(tstrsplit(names(af), '_')[[2]]), 
                      rsum=rsum, asum=asum, tcnt=rsum+asum, af=af)
        af$chr=factor(af$chr, levels=paste0('chr', as.roman(1:16)))

        af%>%filter(tcnt>10) %>% ggplot(aes(x=pos,y=af,color=log2(tcnt))) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
            geom_point() + 
            scale_colour_viridis_c()+
            xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
            facet_grid(~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
            theme_classic()+
            theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
            theme(panel.spacing=unit(0, "lines"))+
            ggtitle(experiment)
        ggsave(paste0(results.dir, 'seg_af.png'), width=22, height=5)
        
        # read in genetic map and format for plotting
        # gmapC=rbindlist(readRDS(paste0(results.dir, 'gmap.RDS')))
        # gmap.subset=gmapC[match(rownames(g.counts[[1]]), paste0(gmapC$chrom, '_', gmapC$ppos)),]
        # gmap.ss=split(gmap.subset, gmap.subset$chrom) 
       
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
      
        # regression based classifier to kick out lousy segs
        uncertainMarkerCount=rowSums(vg<.95 & vg>.05) 
        countsPerCell=colSums(counts)
       
        #png(file=paste0(results.dir, 'additional_seg_filter.png'), width=1024,height=1024)
        plot(log2(countsPerCell), log2(uncertainMarkerCount), 
             xlab='log2(total UMIs per cell)', ylab='log2(uncertain marker count)',main=experiment,sub=ncol(counts)
        )
        seg.classifier.resids=residuals(lm(log2(uncertainMarkerCount)~log2(countsPerCell)))
        outlier.segs=seg.classifier.resids>quantile(seg.classifier.resids, .995)
        ##points(xx[,1][xx3>quantile(xx3,.995)], xx[,2][xx3>quantile(xx3,.995)], col='blue')
        points(log2(countsPerCell)[outlier.segs], log2(uncertainMarkerCount)[outlier.segs], col='blue') 
        #dev.off()

        classifier.name2=names(outlier.segs)[!outlier.segs]
       
        counts=counts[,classifier.name2]
        vg=vg[classifier.name2,]

        cc.big.table=readr::read_delim(state.annot.file, delim='\t')

        cc.df=cc.big.table %>% dplyr::filter( !(cell_cycle %in% cc.bad.vars) & named_dataset == experiment )
      
      #  cc.df$seurat_clusters=as.factor(cc.df$seurat_clusters)
        cc.df=cc.df[cc.df$cell_name %in% rownames(vg),]

        vg=vg[cc.df$cell_name,]
        counts=counts[,cc.df$cell_name]

        exp.results[[as.character(ee)]]$counts=counts
        exp.results[[as.character(ee)]]$vg=vg
        exp.results[[as.character(ee)]]$m.granges=m.granges
        exp.results[[as.character(ee)]]$cell.cycle=cc.df
    }    
    return(exp.results)
}

