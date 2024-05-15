#define a bunch of useful functions and experiment specific variables
source('/data/single_cell_eQTL/yeast/code/processSegregantsSetup.R')

source('/data/single_cell_eQTL/yeast/code/processSegregantsGlobalVars.R')

#load some additional functions for processing the experiment with previously genotyped segregants
source('/data/single_cell_eQTL/yeast/code/processSegregantsPrevGeno.R')

#run HMM and organize data per experiment
source('/data/single_cell_eQTL/yeast/code/processSegregantsGenotyping.R')

#reference 
cell.cycle.assignment.file='/data/single_cell_eQTL/yeast/results/cell_cycle_v5/cell_cycle_feb02162022.tsv'

#combine information from multiple batches for different segregant panels 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    dir.create(comb.out.dir)
    
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

        cc.big.table=readr::read_delim(cell.cycle.assignment.file , delim='\t')
        cc.df=cc.big.table %>% dplyr::filter( cell_cycle != "HALPOID" & cell_cycle != "HAPLOIDS" & named_dataset == experiment )
        cc.df$seurat_clusters=as.factor(cc.df$seurat_clusters)
        cc.df=cc.df[cc.df$cell_name %in% rownames(vg),]

        vg=vg[cc.df$cell_name,]
        counts=counts[,cc.df$cell_name]

        exp.results[[as.character(ee)]]$counts=counts
        exp.results[[as.character(ee)]]$vg=vg
        exp.results[[as.character(ee)]]$m.granges=m.granges
        exp.results[[as.character(ee)]]$cell.cycle=cc.df
    }    

    #combine data sets (investigate size differences (markers) between the different replicates of the 2444 crosses TO DO) 
    vg=do.call('rbind', lapply(exp.results, function(x) x$vg) )
    #rbind(exp.results[[as.character(set.3004[1])]]$vg,exp.results[[as.character(set.3004[2])]]$vg)
    counts=do.call('cbind', lapply(exp.results, function(x) x$counts) )#
    #cbind(exp.results[[as.character(set.3004[1])]]$counts,exp.results[[as.character(set.3004[2])]]$counts)
    m.granges=exp.results[[1]]$m.granges
    cc.df=do.call('rbind', lapply(exp.results, function(x) x$cell.cycle) )
        #rbind(exp.results[[as.character(set.3004[1])]]$cell.cycle,exp.results[[as.character(set.3004[2])]]$cell.cycle)

    #saveRDS(vg, paste0(comb.out.dir, 'vg.RDS'))
    #print(nrow(vg))
    count.non.unique(vg)

    #subset to unique genotypes for cross 3004
    if(set=='3004') {
        bsub=subset.best.unique(vg,counts)
        vg=vg[bsub,]
        counts=counts[,bsub]
        cc.df=cc.df[bsub,]
        count.non.unique(vg)
    }

    pruned=LDprune(vg, m.granges)
    Gsub=pruned$Gsub
    markerGRr=pruned$markerGRr
    count.non.unique(Gsub)


  #  rm(pruned)

    infoCells=rowSums(counts>0)
    expressed.transcripts=names(infoCells)[(infoCells>minInformativeCellsPerTranscript)] #names(sum(infoCells>128)) #[infoCells>128,]

    #select only the transcripts expressed in enough cells 
    Yr=t(counts[expressed.transcripts,])
    transcript.features=data.frame(gene=colnames(Yr), non.zero.cells=colSums(Yr>1), tot.expression=colSums(Yr))

    #cisMarkers=cisMarkers[cisMarkers$transcript %in% expressed.transcripts,]
    #if no subsetting this is not necessary
    #cMsubset=cisMarkers[which(cisMarkers$transcript %in% colnames(Yr)),]

    barcode.features=data.frame(barcode=colnames(counts), nUMI=colSums(counts))

    #replace log(counts)
    mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)+cc.df$named_dataset))
    colnames(mmp1)[3]='experiment'
    
    print("calculating dispersions")

    # custom code for repeated measures mixed models

    if( set == 'Ap') { 
    #if(experiment=="00_BYxRM_480MatA_1" | experiment == "00_BYxRM_480MatA_2") {
       
        segMatch.list=getSegMatch(comb.out.dir, vg, m.granges)
        saveRDS(segMatch.list, file=paste0(comb.out.dir, 'segMatchList.RDS'))

        best_match_seg=segMatch.list$best_match_seg
        gdataPrev=segMatch.list$gdataPrev
        sprevG=segMatch.list$sprevG
        markerGRr2=segMatch.list$markerGRr2
        rm(segMatch.list)

        #run negative binomial mixed model with cell cycle and heritability terms    
        data=data.frame(l_ctot=mmp1[,2], expt=mmp1[,3], Zid=best_match_seg, Cid=cc.df$cell_cycle)
        saveRDS(data, file=paste0(comb.out.dir, 'bGLMM_setup.RDS'))

        #run neg bin mixed model with strain and cell-cycle effects
        dobGLMM(comb.out.dir, data, Yr)

        #additional narrow-sense model for when we have lookup genotypes
        # fit drop1 for each of the random effects in bGLMMs (segregant and cell-cycle), extract VC estimatest, and calc ICCs
        ssH2=calcICCPrevSegs(comb.out.dir, Yr)

        sigY=(which(ssH2[,'H2.q']<.1))

        #fit model with effects of cell-cycle, segregant, and additive effect of all markers
        ssH2a=calcICCPrevSegs_3VC(comb.out.dir ,  Yr,sigY, best_match_seg, sprevG, data) 


            
      df=data.frame(id=factor(best_match_seg))
      mmp2=mmp1[order(df$id),]
      countmat=counts[expressed.transcripts,order(df$id)]
      df.id=df[order(df$id),]
      nullNB=nebula(count=countmat, id=df.id, pred=mmp2) 
    
    } else {
           nullNB=nebula(counts[expressed.transcripts,], mmp1[,1], pred=mmp1)
    }
    

    dispersion.df=data.frame(gene=nullNB$summary$gene, theta.null=1/nullNB$overdispersion[,2])
   
    cisMarkers=getCisMarker(sgd.genes[sgd.genes$gene_id %in% dispersion.df$gene,],
                            markerGRr, t(Yr[,colnames(Yr) %in% dispersion.df$gene]) ) #counts)
    cSplit=split(cisMarkers, cisMarkers$marker)


    if( set == 'Ap') { #experiment=="00_BYxRM_480MatA_1" | experiment == "00_BYxRM_480MatA_2") {
        pGenoCis=doLocalTestPrevGeno(sgd.genes, markerGRr2, Yr, dispersion.df, Gsub, gdataPrev, comb.out.dir, cl)
        clusterExport(cl, varlist=c("mmp1", "Gsub", "counts","best_match_seg","doCisNebula3"))  #, "nebula"))
        system.time({   cisNB=parLapply(cl, cSplit, doCisNebula3)})


    } else {
             clusterExport(cl, varlist=c("mmp1", "Gsub", "counts","doCisNebula2"))  #, "nebula"))
             system.time({   cisNB=parLapply(cl, cSplit, doCisNebula2)})
    }
    names(cisNB)=names(cSplit)
    saveRDS(cisNB, file=paste0(comb.out.dir, 'cisNB2.RDS'))
      
    cis.ps=as.vector(unlist(sapply(cisNB, function(x) x$summary$p_)))
    names(cis.ps)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))

    png(file=paste0(comb.out.dir, 'cis_p_hist.png'), width=512,height=512)
    hist(cis.ps,breaks=50, main=paste(experiment, 'n cells=', nrow(Yr)) , sub=sum(p.adjust(cis.ps,'fdr')<.05))
    dev.off()

    #here we define covariates 
    thetas=as.vector(1/unlist(sapply(cisNB, function(x) x$overdispersion$Cell)))
    names(thetas)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))
    
    cisModel.df=data.frame(gene=names(thetas),theta.cis=thetas)
    dispersion.df=left_join(dispersion.df, cisModel.df, by='gene')
    saveRDS(dispersion.df, file=paste0(comb.out.dir, 'dispersions.RDS'))

    #need to add the marker data structure here 
    segDataList=list(Yr=Yr, 
                        Gsub=Gsub,
                        cisMarkers=cisMarkers, 
                        transcript.features=transcript.features,
                        barcode.features=barcode.features,
                        dispersion.df=dispersion.df,
                        cell.cycle.df=cc.df
    )
    rm(vg)
    rm(exp.results)
    saveRDS(segDataList, file=paste0(comb.out.dir, 'segData.RDS'))

    if(set == 'Ap') {

        #combare local-eQTL betas with lookup and HMM genotypes  -------------------------------
        cisNB=readRDS(paste0(comb.out.dir, 'cisNB2.RDS'))
        cis.ps=as.vector(unlist(sapply(cisNB, function(x) x$summary$p_)))
        cis.beta=as.vector(unlist(sapply(cisNB, function(x) x$summary$logFC_)))
        cis.gene=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))
        nGenoCis=data.frame(gene=cis.gene, beta=cis.beta, p=cis.ps)
        nGenoCis=nGenoCis[nGenoCis$gene %in% pGenoCis$gene,]
    
        geno_cis_comp=list(prevGenoCis=pGenoCis,
             hmmGenoCis=nGenoCis)
        saveRDS(geno_cis_comp, file=paste0(paste0(comb.out.dir, 'geno_cis_comp.RDS')))
        #----------------------------------------------------------------------------------------

       
          }
}


#mapping eQTL within each cell cycle classification for each experiment together 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    barcode.features=segDataList$barcode.features
    dispersion.df=segDataList$dispersion.df

    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene
    
    cc.df=segDataList$cell.cycle.df

    mmp0=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)+cc.df$named_dataset))
    colnames(mmp0)[3]='experiment'

    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cnv=colnames(cc.matrix.manual)
    #  if(set=='3004') { cnv=colnames(cc.matrix.manual)[-1] } else { cnv=colnames(cc.matrix.manual) }
    for(cn in cnv ){
        print(cn)
        cnn=gsub('/',':', cn)
       
        Gsub=segDataList$Gsub[cc.df$cell_name,]
        Yr=segDataList$Yr[cc.df$cell_name,]

        mmp1=mmp0 
        cells.in.cycle=cc.matrix.manual[,cn]==1
        mmp1=mmp1[cells.in.cycle,]
        Gsub=Gsub[cells.in.cycle,]
        Yr=Yr[cells.in.cycle,]

        print('calculating LODs')
        #Yks=Matrix(Yr, sparse=T)
        clusterExport(cl, varlist=c("thetas", "Yr", "Gsub", "mmp1", "nbLL", "domap"))
        clusterEvalQ(cl, { Y=Yr;    DM=mmp1;   return(NULL);})

        bunchoflists=parLapply(cl, names(thetas)[!is.na(thetas)], domap, ss=T)
        LOD=do.call('rbind', lapply(bunchoflists,function(x)x$LRS))
        #LOD=do.call('rbind', parLapply(cl, names(thetas)[!is.na(thetas)], domap) )
        LOD[is.na(LOD)]=0
        LOD=LOD/(2*log(10))
        rownames(LOD)=names(thetas)[!is.na(thetas)]

        saveRDS(LOD, file=paste0(comb.out.dir, 'LOD_NB_', cnn,'.RDS'))

        Betas=do.call('rbind', lapply(bunchoflists,function(x)x$Betas))
        saveRDS(Betas, file=paste0(comb.out.dir, 'Betas_NB_', cnn,'.RDS'))
     
        SEs=do.call('rbind', lapply(bunchoflists,function(x)x$SEs))
        saveRDS(SEs, file=paste0(comb.out.dir, 'SEs_NB_', cnn,'.RDS'))


        print('calculating permutation LODs')
        #permute within batch
        set.seed(100)
        LODpL=list()
        for(i in 1:nperm) { 
            print(i)
            new.index=as.vector(do.call('c', sapply(split(seq_along(mmp1[,'experiment']), as.character(mmp1[,'experiment'])), sample)))
            Ykp=Matrix(Yr[new.index,], sparse=T)
            mmp2=mmp1
            mmp2[,2]=mmp1[new.index,2]
            clusterExport(cl, varlist=c('mmp2', 'Ykp'))
            clusterEvalQ(cl,{ Y=Ykp; DM=mmp2; return(NULL);})
            
            LODp=do.call('rbind', parLapply(cl, names(thetas)[!is.na(thetas)], domap) )
            LODp[is.na(LODp)]=0
            LODp=LODp/(2*log(10))
            rownames(LODp)=names(thetas)[!is.na(thetas)]
            LODpL[[as.character(i)]]=LODp
       }
       LODp=abind(LODpL, along=3)
       rm(LODpL)

       saveRDS(LODp, file=paste0(comb.out.dir, 'LODperm_NB_', cnn,'.RDS')) 
      
       twgFDR= getFDRfx(LOD,LODp, cisMarkers) 
       saveRDS(twgFDR, file=paste0(comb.out.dir, 'fdrfx_NB_', cnn,'.RDS')) #paste0(dout, 'FDRfx.RDS'))

       #markerGR=getMarkerGRanges(list(t(Gsub)))
       #GP=getGlobalPeaks(LOD,markerGR,sgd.genes,fdrfx.fc=twgFDR)
               
       msig=which(LOD>4, arr.ind=T)
       png(file=paste0(comb.out.dir, 'seg_diagnostic_plot_', cnn, '.png'), width=768,height=768)
       plot(msig[,2], msig[,1])
       dev.off()
    }
}
   

#get global peaks and make some plots 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    print(nrow(Yr))
    print(ncol(Yr))
    Gsub=segDataList$Gsub
    print(nrow(Gsub))
    print(ncol(Gsub))


    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    barcode.features=segDataList$barcode.features
    dispersion.df=segDataList$dispersion.df

    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene
    
    cc.df=segDataList$cell.cycle.df
  
    markerGR=getMarkerGRanges(list(t(Gsub)))

    mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)+cc.df$named_dataset))
    colnames(mmp1)[3]='experiment'

   # mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)
    markerGR=getMarkerGRanges(list(t(Gsub)))
   
    L=readRDS(paste0(comb.out.dir,'/LOD_NB_cell_cycle', cycle.cats[1],'.RDS'))
    addL=L
    addL[is.numeric(addL)]=0

    #twgFDR= getFDRfx(LOD,LODp, cisMarkers) 
    #saveRDS(twgFDR, file=paste0(comb.out.dir, 'fdrfx_NB_', cnn,'.RDS')) #paste0(dout, 'FDRfx.RDS'))
    #GP=getGlobalPeaks(LOD,markerGR,sgd.genes,fdrfx.fc=twgFDR)
    
    #fix this 
    addLp=replicate(nperm, addL)   #addLp=abind(addL, addL, addL, addL, addL, along=3)

    BetasList=list()
    SEsList=list()
    cisTableList=list()
    for(cc in cycle.cats) {
        print(cc)
        L=readRDS(paste0(comb.out.dir,'/LOD_NB_cell_cycle', cc,'.RDS'))
        twgFDR=readRDS(file=paste0(comb.out.dir, 'fdrfx_NB_cell_cycle', cc,'.RDS'))

        Betas=readRDS(paste0(comb.out.dir, 'Betas_NB_cell_cycle', cc,'.RDS'))
        SEs=readRDS(paste0(comb.out.dir, 'SEs_NB_cell_cycle', cc,'.RDS'))

        BetasList[[cc]]=Betas
        SEsList[[cc]]=SEs
        GP=getGlobalPeaks(L,markerGR,sgd.genes,fdrfx.fc=twgFDR,Betas=Betas,SEs=SEs)
        saveRDS(GP, file=paste0(comb.out.dir,'/LOD_NB_cell_cycle_peaks_', cc,'.RDS'))

        cmarker.ind=cbind(1:nrow(L),
                    match(cisMarkers$marker[match(rownames(L), cisMarkers$transcript)], colnames(L)) )
        #extract the cis-only scan results 
        cisTableList[[cc]]=data.frame(
            transcript=rownames(L),
            closest.marker=cisMarkers$marker[match(rownames(L), cisMarkers$transcript)],
            LOD=L[cmarker.ind],
            Beta=Betas[cmarker.ind],
            SE=SEs[cmarker.ind] ,
            FDR=twgFDR$c.rtoFDRfx(L[cmarker.ind])
        )

        plotEQTL(GP[GP$FDR<.05,], titleme=cc)
        ggsave(file=paste0(comb.out.dir,'/LOD_NB_', cc, '.png'), width=10, height=10)
      
        Lp=readRDS(paste0(comb.out.dir, 'LODperm_NB_cell_cycle', cc,'.RDS')) 

        addL=addL+L
        addLp=addLp+Lp
    }

   twgFDRG= getFDRfx(addL,addLp, cisMarkers) 
   saveRDS(twgFDRG, file=paste0(comb.out.dir, 'fdrfx_NB_combined.RDS'))
   
   GP=getGlobalPeaks(addL,markerGR,sgd.genes,fdrfx.fc=twgFDRG, Betas=BetasList, SEs=SEsList)
 
   saveRDS(GP, file=paste0(comb.out.dir,'/LOD_NB_combined_peaks.RDS'))
   plotEQTL(GP[GP$FDR<.05,], titleme='combined', CI=F)
   ggsave(file=paste0(comb.out.dir,'/LOD_NB_combined.png'), width=10, height=10)
   
   cmarker.ind=cbind(1:nrow(addL),
                    match(cisMarkers$marker[match(rownames(addL), cisMarkers$transcript)], colnames(addL)) )

    cisTableList[['combined']]=data.frame(
            transcript=rownames(addL),
            closest.marker=cisMarkers$marker[match(rownames(addL), cisMarkers$transcript)],
            LOD=addL[cmarker.ind],
            #Beta=Betas[cmarker.ind],
            #SEs=SEs[cmarker.ind] ,
            FDR=twgFDRG$c.rtoFDRfx(addL[cmarker.ind])
        )
     saveRDS(cisTableList, file=paste0(comb.out.dir,'/cis_only_test_CombinedResults.RDS'))

  str(cisTableList) 
   print(nrow(GP[GP$FDR<.1,]))
   print(nrow(GP[GP$FDR<.05,]))
}

# cell cycle x cis interactions 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')

    cisTableList=readRDS(file=paste0(comb.out.dir,'/cis_only_test_CombinedResults.RDS'))

    # test for sig interactions between cell cycle betas
    #process cis_only_test_CombinedResults.RDS
    cTL=cisTableList[-length(cisTableList)]
    ucontrasts=combn(names(cTL),2)
    cTL.contrasts=list()
    for (i in 1:ncol(ucontrasts)) {
        coi=paste(ucontrasts[1,i], ucontrasts[2,i])

        B1=cTL[[ucontrasts[1,i]]]$Beta
        B2=cTL[[ucontrasts[2,i]]]$Beta

        SE1=cTL[[ucontrasts[1,i]]]$SE
        SE2=cTL[[ucontrasts[2,i]]]$SE

        Zstat=(B1-B2)/sqrt(SE1^2+SE2^2)
        pval=2*pnorm(-abs(Zstat),lower.tail=T)
        cTL.contrasts[[coi]]=data.frame(transcript=cTL[[ucontrasts[1,i]]]$transcript, Zstat=Zstat, pval=pval)
    }
    saveRDS(cTL.contrasts, file=paste0(comb.out.dir,'/cis_only_test_CombinedResultsContrasts.RDS'))

    sum(p.adjust(unlist(lapply(cTL.contrasts, function(x) x$pval)), method='fdr')<.05)
    sort(unique(unlist(lapply(cTL.contrasts, function(x) x$transcript))[which(p.adjust(unlist(lapply(cTL.contrasts, function(x) x$pval)), method='fdr')<.05)]))


}



















# summary of QTL detected at different FDR thresholds 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')

    GP=readRDS(paste0(comb.out.dir,'/LOD_NB_combined_peaks.RDS'))
    print('5% FDR')
    print(nrow(GP[GP$chrom==GP$tchr & GP$FDR<.05,]))
    print(nrow(GP[GP$chrom!=GP$tchr & GP$FDR<.05,]))

    print('10% FDR')
    print(nrow(GP[GP$chrom==GP$tchr & GP$FDR<.1,]))
    print(nrow(GP[GP$chrom!=GP$tchr & GP$FDR<.1,]))


}

# mapping cell cycle classifications 
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    barcode.features=segDataList$barcode.features
    dispersion.df=segDataList$dispersion.df

    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene
    
    cc.df=segDataList$cell.cycle.df

    mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)+cc.df$named_dataset))
    colnames(mmp1)[3]='experiment'

   # mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)

    clusterExport(cl, varlist=c("Gsub", "mmp1", "cc.incidence", "domap_logistic"))
    #clusterEvalQ(cl, { Y=Yr;    DM=mmp1;   return(NULL);})
    bunchoflists=parLapply(cl, colnames(cc.incidence), domap_logistic, ss=T)

    LOD=do.call('rbind', lapply(bunchoflists, function(x)x$LRS))
    LOD[is.na(LOD)]=0
    LOD=LOD/(2*log(10))
    rownames(LOD)=colnames(cc.incidence)
    #dir.create(paste0(base.dir, 'results/cell_cycle_v5/', experiment,'/'))

    Gsub.s=scale(Gsub)
    Cmat=(t(Gsub.s) %*% Gsub.s)/(nrow(Gsub)-1)
    meff=getMeff_Li_and_Ji(Cmat)
    #attr(LOD,'FWER_thresh_1%')= qchisq(2*(.01/meff),1,lower.tail=F)/(2*log(10))
    attr(LOD,'FWER_thresh_5%')= qchisq(2*(.05/(meff*5)),1,lower.tail=F)/(2*log(10))
    print(attr(LOD, 'FWER_thresh_5%'))
    saveRDS(LOD, file=paste0(comb.out.dir, 'cell_cycle_assignment_LOD.RDS'))
        
    #qchisq(2*(.05/meff),1,lower.tail=F)/(2*log(10))

    saveRDS(LOD, file=paste0(comb.out.dir, 'cell_cycle_assignment_LOD.RDS'))

    Betas=do.call('rbind', lapply(bunchoflists,function(x)x$Betas))
    rownames(Betas)=rownames(LOD) #colnames(cc.incidence)

    saveRDS(Betas, file=paste0(comb.out.dir, 'cell_cycle_assignment_Betas.RDS'))
     
    SEs=do.call('rbind', lapply(bunchoflists,function(x)x$SEs))
    rownames(SEs)=rownames(LOD)
    saveRDS(SEs, file=paste0(comb.out.dir, 'cell_cycle_assignment_SEs.RDS'))

    pdf(file=paste0(comb.out.dir, 'cell_cycle_assignment_LOD.pdf'), width=10, height=5)
    for(i in 1:nrow(LOD)){
        plot(LOD[i,],main=rownames(LOD)[i], ylab='LOD', xlab='marker index')
        abline(v=cumsum(c(0,rle(tstrsplit(colnames(Gsub), '_')[[1]])$lengths)), lty=2, col='blue')
    }
    dev.off()
}


   #read cell cycle LODs
   
for(set in names(sets)){
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')

    ccLOD=readRDS(file=paste0(comb.out.dir, 'cell_cycle_assignment_LOD.RDS'))
    ccLOD=ccLOD[1:5,]
    rownames(ccLOD)=gsub('cell_cycle', '', rownames(ccLOD))
    #tidyr::gather(t(data.frame(ccLOD)), key="cell_cycle", value="LOD")
    df=data.frame(chr=as.character(seqnames(markerGR)), pos=start(markerGR), t(ccLOD))
    df=tidyr::gather(df, key="cell_cycle", value="LOD", -chr,-pos)


    ccLOD=ggplot(df,aes(pos,LOD)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
        geom_line()+
        xlab('')+ylab('')+
        # scale_alpha(guide = 'none') + 
        facet_grid(cell_cycle~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        #theme(panel.spacing=unit(0, "lines"))+
        geom_hline(aes(yintercept=4, colour='red'))+theme(legend.position='none')+
        ggtitle(set)
   plot(ccLOD) 
   ggsave(file=paste0(comb.out.dir,'/cell_cycle_assignment_LOD.png'), width=10, height=8)
    
}

#get significant hotspot bins and test for trans-eQTL x CC interactions 
source('/data/single_cell_eQTL/yeast/code/processSegregantsHotspots.R')



#}

#epistasis -----------------------------------------------------------------------
doInt=function(gn, ...) { 
    print(gn)
    theta.est=thetas[gn]
    YY=Yr[,gn]
    nbn=negative.binomial(theta=theta.est,link='log')
          
    GPsdf=GPs[[gn]]
    #Xmat.0=cbind(mmp2, cc.matrix.manual)
    peak.combos=expand.grid(GPsdf$peak.marker,GPsdf$peak.marker)
    peak.combos=peak.combos[peak.combos[,1]!=peak.combos[,2],]
    peak.combos[,1]=as.character(peak.combos[,1])
    peak.combos[,2]=as.character(peak.combos[,2])

    LRS2=rep(NA, nrow(peak.combos))

    for(gx in 1:nrow(peak.combos)){
      print(gx)
      DM=cbind(mmp2, Gsub[,peak.combos[gx,1]],  Gsub[,peak.combos[gx,2]])
      fnbrN=fastglmPure(DM, YY, family=nbn)
      nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)
      DMint=cbind(DM,  Gsub[,peak.combos[gx,1]]* Gsub[,peak.combos[gx,2]])
      fnbrF=fastglmPure(DMint, YY,  family=nbn)
      LRS2[gx]=-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))
    }
    results=(cbind(peak.combos, LRS2))
    colnames(results)=c('Q1', 'Q2', 'LRS')
    return(results)
}

#epistasis -----------------------------------------------------------------------
for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    barcode.features=segDataList$barcode.features
    dispersion.df=segDataList$dispersion.df

    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene
    
    cc.df=segDataList$cell.cycle.df

    mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)+cc.df$named_dataset))
    colnames(mmp1)[3]='experiment'

   # mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)
    markerGR=getMarkerGRanges(list(t(Gsub)))
    GP =  readRDS(paste0(comb.out.dir,'/LOD_NB_combined_peaks.RDS'))

    GPs=split(GP, GP$transcript)
   # gns=names(GPs)
   # gn=gns[1]
    
    cnv=colnames(cc.matrix.manual)
    #  if(set=='3004') { cnv=colnames(cc.matrix.manual)[-1] } else { cnv=colnames(cc.matrix.manual) }
    for(cn in cnv ){
        #this 
        cnn=gsub('/',':', cn)

        print(cn)
        Gsub=segDataList$Gsub
        Yr=segDataList$Yr

        mmp2=mmp1 
        cells.in.cycle=cc.matrix.manual[,cn]==1
        mmp2=mmp2[cells.in.cycle,]
        Gsub=Gsub[cells.in.cycle,]
        Yr=Yr[cells.in.cycle,]

        print('calculating QTLxQTL LODs')
        #Yks=Matrix(Yr, sparse=T)
        clusterExport(cl, varlist=c("thetas", "GPs", "Yr", "Gsub", "mmp2", "nbLL", "doInt"))
        clusterEvalQ(cl, { Y=Yr;    DM=mmp2;   return(NULL);})


       #perform 2D scan between chr max ------------------------------------
       QQ =   parLapply(cl, names(GPs)[1:16], doInt) 
       names(QQ)=names(GPs)[1:16]
       QQ=data.table::rbindlist(QQ, idcol='transcript')
       QQ$LOD=QQ$LRS
       QQ$LOD[is.na(QQ$LOD)]=0
       QQ$LOD=QQ$LOD/(2*log(10))
       saveRDS(QQ, file=paste0(comb.out.dir, 'LOD_NB_QQ_', cnn,'.RDS'))
       #--------------------------------------------------------------------

   #permutations        ------------------------------------------------------------------------------------------------    
        set.seed(100)
        LODpL=list()
        for(i in 1:nperm) { 
            print(i)
            new.index=as.vector(do.call('c', sapply(split(seq_along(mmp2[,'experiment']), as.character(mmp2[,'experiment'])), sample)))
            Ykp=Matrix(Yr[new.index,], sparse=T)
            mmp3=mmp2
            mmp3[,2]=mmp2[new.index,2]
            clusterExport(cl, varlist=c('mmp3', 'Ykp'))
            clusterEvalQ(cl,{ Yr=Ykp; DM=mmp3; return(NULL);})
            
            QQp =   parLapply(cl, names(GPs), doInt) 
            names(QQp)=names(GPs) 
            QQp=data.table::rbindlist(QQp, idcol='transcript')
            QQp$LOD=QQp$LRS
            QQp$LOD[is.na(QQp$LOD)]=0
            QQp$LOD=QQp$LOD/(2*log(10))
            LODpL[[i]]=QQp
        }

        LODpLn=LODpL[[1]][,-4]
        for( i in 2:nperm) {  LODpLn=cbind(LODpLn, LODpL[[i]]$LOD) }
        saveRDS(LODpLn, file=paste0(comb.out.dir, 'LODperm_NB_QQ_', cnn, '.RDS'))
  #---------------------------------------------------------------------------------------------------------------------
    
    }
}


#process epistasis

iLODs=list()

for(set in names(sets)){
    print(set)
    comb.out.dir=paste0(base.dir, 'results/combined/', set, '/')
    segDataList=readRDS(paste0(comb.out.dir, 'segData.RDS'))
    cc.df=segDataList$cell.cycle.df

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    barcode.features=segDataList$barcode.features
    dispersion.df=segDataList$dispersion.df

    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene
 
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)
    markerGR=getMarkerGRanges(list(t(Gsub)))

    cnv=colnames(cc.matrix.manual)
    #  if(set=='3004') { cnv=colnames(cc.matrix.manual)[-1] } else { cnv=colnames(cc.matrix.manual) }
    QQs=list()
    QQsp=list()
    for(cn in cnv ){
        print(cn)
        cnn=gsub('/',':', cn)
        QQtemp=readRDS(file=paste0(comb.out.dir, 'LOD_NB_QQ_', cnn,'.RDS'))
        Q1ind=match(QQtemp$Q1, markerGR$name)
        Q2ind=match(QQtemp$Q2, markerGR$name)

        Q1f=ifelse(Q1ind<Q2ind, Q1ind,  Q2ind)
        Q2f=ifelse(Q2ind>Q1ind, Q2ind,  Q1ind)
        QQtemp$Q1= markerGR$name[Q1f]
        QQtemp$Q2= markerGR$name[Q2f]
        ddQ=paste0(QQtemp$transcript, QQtemp$Q1, QQtemp$Q2)
        QQs[[cnn]]=QQtemp[!duplicated(ddQ),]
        QQtemp.perm=readRDS(paste0(comb.out.dir, 'LODperm_NB_QQ_', cnn, '.RDS'))
        QQsp[[cnn]]=QQtemp.perm[!duplicated(ddQ),]
    }
    iLODs[[set]]$obs=QQs
    iLODs[[set]]$exp=QQsp
}


lapply(iLODs, function(x) lapply(x, function(y) y[y$LOD>6,]))


#process stats for observed data 
jLOD=lapply(iLODs, function(QQs) {

    jLOD=rowSums(sapply(QQs$obs, function(x) x$LOD))
    jLOD=data.frame(QQs$obs[[1]][,c(1:3)], LOD=jLOD)
    return(jLOD)                             } )

#process stats for perm data
jLODp=lapply(iLODs, function(QQs) {

    jLODt=rowSums(do.call('cbind', (lapply(QQs$exp, function(x) x[,4]))))
    jLODt=cbind(jLODt, rowSums(do.call('cbind', (lapply(QQs$exp, function(x) x[,5])))))
    jLODt=cbind(jLODt, rowSums(do.call('cbind', (lapply(QQs$exp, function(x) x[,6])))))
    jLODt=cbind(jLODt, rowSums(do.call('cbind', (lapply(QQs$exp, function(x) x[,7])))))
    jLODt=cbind(jLODt, rowSums(do.call('cbind', (lapply(QQs$exp, function(x) x[,8])))))
   
 
    jLOD=data.frame(QQs$exp[[1]][,c(1:3)], jLODt) 
    
    #LOD=jLOD)
    return(jLOD)                             } )


#let's just do FDR manually here 

for(set in names(jLOD) { 
   #    set

    x=jLOD[[set]]
    max.obsLOD=sapply(split(x$LOD, x$transcript), max)

    vint=seq(0,max(max.obsLOD)+.01, .01)

    x=jLODp[[set]]
    exp=sapply(4:8, function(i) {
      sapply(split(x[,i], x$transcript), max) })

    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = vint

    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return((apply(exp,2,function(x) sum(x>thresh))))
                 })
    expPcnt= apply(expPcnt, 2, mean)
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    
    pFDR = rev(cummax(rev(pFDR)))
    g.fdrFX=approxfun(pFDR, vint, ties=mean, rule=2)
    g.rtoFDRfx=approxfun(vint,pFDR, ties=mean, rule=2)

    }



    
lapply(jLOD, function(x) x[x$LOD>6,])

hist(pchisq((2*log(10))*jLOD[[3]]$LOD,1, lower.tail=F))




