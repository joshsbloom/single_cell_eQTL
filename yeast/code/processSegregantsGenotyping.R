#build sets of recurring hets and distorted/unreliable markers  
for(setn in names(sets)){

    print(setn)
    
    exp.results=list()

    # iterate over each experiment -----------------------------------------------
    for( ee in  sets[[setn]] ) { 

    experiment=experiments[ee]
    cross=crossL[[cdesig[ee]]]
    parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')
    dir.create(results.dir)

    #clean up this crap-------------------------
    gmap=genetic.maps[[cdesig[ee]]]
    gmap.s=gmap; colnames(gmap.s)[2]='map'
    gmap.s=split(gmap.s, gmap.s$chrom)
    gmap.s=gmap.s[chroms]
    gmap.s=jitterGmap(gmap.s)
    saveRDS(gmap.s, file=paste0(results.dir, 'gmap.RDS'))
    #-------------------------------------------
   
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
        saveRDS(counts, file=paste0(results.dir, 'counts.RDS'))
    #-------------------------------------------------------------

    # counts per parental allele (ref.counts and alt.counts, names not accurate) (fixes phasing) --------
    # also filters out counts from heavily distorted variants /// update this with diploid ASE data or parental data when we have it
        g.counts=getGenoInformativeCounts(data.dir, cranger.dir, gmap,parents, log2diff=8)
        saveRDS(g.counts, file=paste0(results.dir, 'gcounts.RDS'))
    #-----------------------------------------------------------------------------------------------------
    
    # can filter here for only sites that have informative variants in this data set,
    # but this may make comparison across datasets difficult in the future, hold on this for now
        rQTL.coded=encodeForRQTL(g.counts)
        #saveRDS(rQTL.coded, file=paste0(results.dir, 'rQTLcoded.RDS'))
    #----------------------------------------------------------------------------------------------------

    het.count=apply(rQTL.coded, 2, function(x) sum(x==0,na.rm=T))
    tot.count=colSums(g.counts[[1]]+g.counts[[2]]) #ref.counts+alt.counts)
   
    plot(log2(tot.count),log2(het.count),col="#00000022")  # manual inspection of total counts vs sites classified as hets
    
    # note this is classification as probably haploid (TRUE) or not (dipoid, remated segregant/doublet FALSE)
    classification=het.count<het.thresholds[ee]
    saveRDS(classification, file=paste0(results.dir, 'cellFilter.RDS'))

    png(file=paste0(results.dir, 'plot_classify_haploids.png'), width=1024, height=1024)
    ## add information about number of cells classified as haploid vs not
    plot(log2(tot.count),log2(het.count), col=  classification+1,
         xlab='log2(total geno informative reads)', ylab='log2(total heterozygous sites)', 
         main=experiment, sub=paste('total = ', length(het.count),  ' haploid =', sum(classification) ) )#s[[experiment]])))
    dev.off()

    # identify unreliable sites (sites called as hets in some fraction of cells, re-evaluate threshold chosen here)
    het.cells=!classification
    recurring.hets=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==0,na.rm=T)) > (het.across.cells.threshold*sum(!het.cells))
    print(sum(recurring.hets))
    print(sum(recurring.hets)/nrow(rQTL.coded)) 

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
    #af[(rsum+asum)<10]=NA
    af=data.frame(chr=tstrsplit(names(af), '_')[[1]], pos=as.numeric(tstrsplit(names(af), '_')[[2]]), 
                     rsum=rsum, asum=asum, tcnt=rsum+asum, af=af)
    af$chr=factor(af$chr, levels=paste0('chr', as.roman(1:16)))
    af$af.folded=abs(af$af-.5)
    af$recurring.hets=recurring.hets
    rownames(af)=names(recurring.hets)

    saveRDS(af, file=paste0(results.dir, 'af.RDS'))
    rm(rQTL.coded)

    }
}




    #per suggestion of James, run HMM on everything to start
    
    #cross2=buildCross2(rQTL.coded, gmap, 'riself', het.sites.remove=T, het.cells.remove=F,
    #                   het.cells=het.cells, recurring.hets=recurring.hets)
    #rqtl.genoprobs=calc_genoprob(cross2, error_prob=.005, lowmem=F, quiet=F, cores=16)
    #rqtl.genoprobs=do.call('cbind', sapply(rqtl.genoprobs, function(x) x[,2,]))
    #saveRDS(rqtl.genoprobs, file=paste0(results.dir, 'rqtl_genoprobs.RDS'))
    #emissionProbs=estimateEmissionProbs(g.counts[[1]],  g.counts[[2]], error.rate=.005,
    #                                    recurring.het.sites.remove=T,
    #                                    recurring.hets=recurring.hets)
    #test=runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms[3], calc='viterbi',return.chr=T) #, n.indiv=1000)
    #par(yaxs='i')
    #plot(colSums(test[which(classification),]-1)/sum(classification), ylim=c(0,1))
    #abline(h=.85)

#aggregate results 
af.results=list()
for(setn in names(sets)){

    print(setn)
    comb.out.dir=paste0(base.dir, 'results/combined/', setn, '/')
    dir.create(comb.out.dir)
    

    # iterate over each experiment -----------------------------------------------
    for( ee in  sets[[setn]] ) { 

    #for(ee in c(1:6,13:16) ) { #5:length(experiments)){
    experiment=experiments[ee]
    cross=crossL[[cdesig[ee]]]
    parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')
     af.results[[setn]][[experiment]]=readRDS(file=paste0(results.dir, 'af.RDS'))
     af.results[[setn]][[experiment]]$flagged.marker=af.results[[setn]][[experiment]]$af.folded>.4 & af.results[[setn]][[experiment]]$tcnt>10
    }
}


#assemble bad markers 
bad.marker.list=lapply(af.results, function(y) { 
             b=rowSums(sapply(y, function(x) { return(x$recurring.hets | x$flagged.marker) } ) )
             names(b)=rownames(y[[1]])
             return(b)
                     })
#saveRDS(bad.marker.list, file=paste0(base.dir, 'results/badMarkerList.RDS'))

#calculate hmm
for(setn in names(sets)){

    print(setn)
    comb.out.dir=paste0(base.dir, 'results/combined/', setn, '/')
    dir.create(comb.out.dir)
    
    exp.results=list()

    recurring.hets=bad.marker.list[[setn]]
    recurring.hets=recurring.hets>0

    # iterate over each experiment -----------------------------------------------
    for( ee in  sets[[setn]] ) { 
        print(ee)
        #for(ee in c(1:6,13:16) ) { #5:length(experiments))
        experiment=experiments[ee]
        cross=crossL[[cdesig[ee]]]
        parents=parentsL[[cdesig[ee]]]
        data.dir=data.dirs[ee]
        results.dir=paste0(base.dir, 'results/', experiment, '/')


        gmap.s=readRDS(paste0(results.dir, 'gmap.RDS'))
        
        g.counts=readRDS(paste0(results.dir, 'gcounts.RDS'))
        tmp=as_matrix(g.counts$ref.counts[,-recurring.hets])
        #rsum=Rfast::colsums(tmp)
        #rm(tmp)
        tmp1=as_matrix(g.counts$alt.counts[,-recurring.hets])
        asum=Rfast::rowsums(tmp)
        #rm(tmp1)
        tmp2=tmp+tmp1
        print(colsums(tmp2>0))
        #note recurring.hets contains all flagged markers for removal
        emissionProbs=estimateEmissionProbs(g.counts[[1]],  g.counts[[2]], error.rate=.005,
                                            recurring.het.sites.remove=T,
                                            recurring.hets=recurring.hets)
        print('calculating genotype probabilities')
        runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='genoprob') #, n.indiv=1000)
        print('calculating viterbi path')
        runHMM(emissionProbs[[1]], emissionProbs[[2]], results.dir,gmap.s, chroms, calc='viterbi') #, n.indiv=1000)
}

}
 
