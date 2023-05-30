#define a bunch of useful functions and experiment specific variables
source('/data/single_cell_eQTL/yeast/code/processSegregantsSetup.R')

#load some additional functions for processing the experiment with previously genotyped segregants
source('/data/single_cell_eQTL/yeast/code/processSegregantsPrevGeno.R')

#run HMM and organize data per experiment
source('/data/single_cell_eQTL/yeast/code/processSegregantsGenotyping.R')


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

        cc.big.table=readr::read_delim('/data/single_cell_eQTL/yeast/results/cell_cycle_v5/cell_cycle_feb02162022.tsv', delim='\t')
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

    pruned=LDprune(vg, m.granges)
    Gsub=pruned$Gsub
    markerGRr=pruned$markerGRr
  
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








#not incorporated yet 03/29/22
fitNBTransModel=function(tcP, Yk, Gsub,  cisNB, mmp1, cMsubset,sig=.1 ) {
        tcPsig=data.frame(tcP[tcP$FDR<sig,])
        tcPsig$negbin.beta=NA
        tcPsig$negbin.se=NA
        #tcPsig$LLik=NA
        #tcPsig$negbin.LRS=NA
        tcPsig$negbin.p=NA
        #mmp0=mmp1[,-1]
        for(r in 1:nrow(tcPsig)){
            print(r)
            gn=as.character(tcPsig[r,]$transcript)
            #pcis=as.character(cMsubset$marker[match(gn, cMsubset$transcript)])
            pmarker=as.character(tcPsig[r,]$peak.marker)
            theta.est=cisNB[[gn]]$theta
            #nmLLik=cisNB[[gn]]$LLik
            XX=cbind(mmp1, Gsub[,pmarker])
            msize=ncol(XX)
            fnbrF=fastglm(XX, Yk[,gn],family=negative.binomial(theta=theta.est,link='log'), maxit=500)
            tcPsig$negbin.beta[r]=as.numeric(coef(fnbrF)[msize])
            tcPsig$negbin.se[r]=as.numeric(fnbrF$se[msize])
            tcPsig$negbin.p[r]=pchisq((2*log(10))*tcPsig[r,]$LOD,1, lower.tail=F)
            #tcPsig$LLik[r]=as.numeric(logLik(fnbrF))
            #tcPsig$negbin.LRS[r]=as.numeric(-2*(nmLLik-logLik(fnbrF)))
            #pchisq( tcPsig$negbin.LRS[r],1, lower.tail=F) #-2*(nmodel-plr)
            print(tcPsig[r,])
        }
        return(tcPsig)
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



















#    plotEQTL(GP[GP$LOD>4,], titleme='combined', CI=F)
#    ggsave(file=paste0(comb.out.dir,'/LOD_NB_combined.png'), width=10, height=10)
#

#sgd.genes=sgd.granges[sgd.granges$type=='gene',]



#ByxRM previously genotypoed




#for(ee in c(2,4,7:12,15,16)) { #1:length(experiments)){}}
#for(ee in c(15,16)) { #1:length(experiments)){}}
#c(1:4,7:10,11,12)) { #,15,16)) { #1:length(experiments)){}}
for(ee in c(13,14) ) { 

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
    png(file=paste0(results.dir, 'additional_seg_filter.png'), width=1024,height=1024)
    plot(log2(countsPerCell), log2(uncertainMarkerCount), 
         xlab='log2(total UMIs per cell)', ylab='log2(uncertain marker count)',main=experiment,sub=ncol(counts)
    )
    seg.classifier.resids=residuals(lm(log2(uncertainMarkerCount)~log2(countsPerCell)))
    outlier.segs=seg.classifier.resids>quantile(seg.classifier.resids, .995)
    #points(xx[,1][xx3>quantile(xx3,.995)], xx[,2][xx3>quantile(xx3,.995)], col='blue')
    points(log2(countsPerCell)[outlier.segs], log2(uncertainMarkerCount)[outlier.segs], col='blue') 
    dev.off()

    classifier.name2=names(outlier.segs)[!outlier.segs]
    counts=counts[,classifier.name2]
    vg=vg[classifier.name2,]

    #code to match original BYxRM with previous genotypes -------------------
    # load previous data 
    if(experiment=="00_BYxRM_480MatA_1" | experiment == "00_BYxRM_480MatA_2") {
        load('/data/eQTL/gdata_42k.RData')
        prevG=gdata
        rm(gdata)
        colnames(prevG)=gsub(':', '_', colnames(prevG))
        scname=tstrsplit(colnames(prevG),'_')
        scname=paste0(scname[[1]], '_', scname[[2]], '_', tstrsplit(scname[[3]],'/')[[1]], '_',tstrsplit(scname[[3]],'/')[[2]] )
        colnames(prevG)=scname

    # line up with existing genotypes 
        vg2=vg
        vg2name=tstrsplit(colnames(vg2), '_')
        colnames(vg2)=paste0(vg2name[[1]], '_', vg2name[[2]], '_', vg2name[[3]], '_', vg2name[[4]])

        matchingmarkers=colnames(vg2)[colnames(vg2) %in% colnames(prevG)]
        vg2=vg2[,matchingmarkers]
        prevG=prevG[,matchingmarkers]
        stvg2=scale(t(vg2))
        sprevG=scale(t(prevG))   
    # find best match genotype
        allcors=crossprod(stvg2, sprevG)/(nrow(sprevG)-1)
        best_match_seg=rownames(prevG)[(apply(allcors, 1 , which.max))]
        best_match_seg.r=apply(allcors, 1 , max)
        rm(prevG); rm(vg2); rm(stvg2); rm(sprevG)
    }

    #plot(log2(rowSums(Yr)), uncertainMarkerCount)
    #xx=cbind(log2(rowSums(Yr)), log2(uncertainMarkerCount))
    #plot(xx[,1], xx[,2], xlab='log2(total UMIs per cell)', ylab='log2(uncertain marker count))
    #xx3=residuals(lm(xx[,2]~xx[,1]))
    #points(xx[,1][xx3>quantile(xx3,.995)], xx[,2][xx3>quantile(xx3,.995)], col='blue')

    pruned=LDprune(vg, m.granges)
    Gsub=pruned$Gsub
    markerGRr=pruned$markerGRr
    rm(pruned)

    # af diagnostic    
    print("making diagonstic plots")
    png(file=paste0(results.dir, 'seg_diagnostic_af_plot.png'), width=1024,height=1024)
    par(mfrow=c(2,1))
    plot(m.granges$gcoord, colSums(vg)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='fraction of segs that are parent 2', 
         main=experiment, sub=paste(ncol(vg), 'markers'))
    abline(v=gcoord.key, lty=2, col='lightblue')
    abline(h=0.5, col='grey')
    #dev.off()
    #png(file=paste0(results.dir, 'seg_diagnostic_af_plot_pruned.png'), width=1024,height=512)
    plot(markerGRr$gcoord, colSums(Gsub)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='fraction of segs that are parent 2',
         main=experiment, sub=paste(ncol(Gsub), 'markers after pruning'))
    abline(v=gcoord.key, lty=2, col='lightblue')
    abline(h=0.5, col='grey')
    dev.off()
   
    rm(vg) 
    # relatedness between segs 
    tG=t(Gsub)
    tt=hclust(as.dist(1-(crossprod(scale(t(Gsub)))/(nrow(tG)-1))^2))
    cut.tt=cutree(tt,h=.5)
    cut.tt.table=table(cut.tt)
    singleton.cells=names(cut.tt)[cut.tt %in% as.numeric(names(cut.tt.table[cut.tt.table==1]))]
    tt2=hclust(as.dist(1-(crossprod(scale(t(Gsub[singleton.cells,])))/(nrow(tG)-1))^2))
    plot(tt2, labels=F, main=experiment, ylab="1-r^2")

    png(file=paste0(results.dir, 'seg_diagnostic_corr_clust_plot.png'), width=1024,height=512)
    plot(tt, labels=F, main=experiment, ylab="1-r^2")
    dev.off()
    
    #count informative cells per transcript
    infoCells=rowSums(counts>0)
    expressed.transcripts=names(infoCells)[(infoCells>minInformativeCellsPerTranscript)] #names(sum(infoCells>128)) #[infoCells>128,]

    #select only the transcripts expressed in enough cells 
    Yr=t(counts[expressed.transcripts,])
    
    cisMarkers=getCisMarker(sgd.genes[sgd.genes$gene_id %in% rownames(counts),], markerGRr, t(Yr)) #counts)
    #cisMarkers=cisMarkers[cisMarkers$transcript %in% expressed.transcripts,]
    #if no subsetting this is not necessary
    #cMsubset=cisMarkers[which(cisMarkers$transcript %in% colnames(Yr)),]

    transcript.features=data.frame(gene=colnames(Yr), non.zero.cells=colSums(Yr>1), tot.expression=colSums(Yr))
    barcode.features=data.frame(barcode=colnames(counts), nUMI=colSums(counts))

    #replace log(counts)
    mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)))
    cSplit=split(cisMarkers, cisMarkers$marker)
   
    print("calculating dispersions")
    if(experiment=="00_BYxRM_480MatA_1" | experiment == "00_BYxRM_480MatA_2") {
         clusterExport(cl, varlist=c("mmp1", "Gsub", "counts","best_match_seg","doCisNebula3"))  #, "nebula"))
         system.time({   cisNB=parLapply(cl, cSplit, doCisNebula3)})

    } else {
         clusterExport(cl, varlist=c("mmp1", "Gsub", "counts","doCisNebula2"))  #, "nebula"))
         system.time({   cisNB=parLapply(cl, cSplit, doCisNebula2)})
   
         #982
    }
    names(cisNB)=names(cSplit)


    saveRDS(cisNB, file=paste0(results.dir, 'cisNB2.RDS'))
    cis.ps=as.vector(unlist(sapply(cisNB, function(x) x$summary$p_)))
    names(cis.ps)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))

    png(file=paste0(results.dir, 'cis_p_hist.png'), width=512,height=512)
    hist(cis.ps,breaks=50, main=paste(experiment, 'n cells=', nrow(Yr)) , sub=sum(p.adjust(cis.ps,'fdr')<.05))
    dev.off()

    #here we define covariates 
    thetas=as.vector(1/unlist(sapply(cisNB, function(x) x$overdispersion$Cell)))
    names(thetas)=as.vector(unlist(sapply(cisNB, function(x) x$summary$gene)))
    
    if(experiment=="00_BYxRM_480MatA_1" | experiment == "00_BYxRM_480MatA_2") {
          
        df=data.frame(id=factor(best_match_seg))
        mmp2=mmp1[order(df$id),]
        countmat=counts[expressed.transcripts,order(df$id)]
        df.id=df[order(df$id),]
        nullNB=nebula(count=countmat, id=df.id, pred=mmp2) 
    } else {
        nullNB=nebula(counts[expressed.transcripts,], mmp1[,1], pred=mmp1)
    }

    dispersion.df=data.frame(gene=nullNB$summary$gene, theta.null=1/nullNB$overdispersion[,2])
    cisModel.df=data.frame(gene=names(thetas),theta.cis=thetas)
    dispersion.df=left_join(dispersion.df, cisModel.df, by='gene')
    saveRDS(dispersion.df, file=paste0(results.dir, 'dispersions.RDS'))

    #need to add the marker data structure here 
    segDataList=list(Yr=Yr, 
                    Gsub=Gsub,
                    cisMarkers=cisMarkers, 
                    transcript.features=transcript.features,
                    barcode.features=barcode.features,
                    dispersion.df=dispersion.df)

    saveRDS(segDataList, file=paste0(results.dir, 'segData.RDS'))
}












































# mapping within each cell cycle category
for (ee in c(1:4,7:12) ) {
    experiment=experiments[ee]
    print(ee)
    print(experiment)    
    #cross=crossL[[cdesig[ee]]]
    #parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')

    segDataList=readRDS(file=paste0(results.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    #barcode.featues=segDataList$barcode.features

    dispersion.df=segDataList$dispersion.df
    counts=readRDS(paste0(results.dir, 'counts.RDS'))

    # replace if barcode.features is present
    #mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)))
  
    #could switch between thetas here 
    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene


    cc.big.table=readr::read_delim('/data/single_cell_eQTL/yeast/results/cell_cycle_v5/cell_cycle_feb02162022.tsv', delim='\t')
    cc.df=cc.big.table %>% dplyr::filter( cell_cycle != "HALPOID" & cell_cycle != "HAPLOIDS" & named_dataset == experiment )
    #cell_cycle != "ALPHA" &
   
    cc.df=cc.df[cc.df$cell_name %in% rownames(Gsub),]
    cc.df$seurat_clusters=as.factor(cc.df$seurat_clusters)
   
    Gsub=Gsub[cc.df$cell_name,]
    Yr=Yr[cc.df$cell_name,]

    print(all.equal(cc.df$cell_name, rownames(Gsub)))

    mmp0=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    
    #model per cell cycle classification
    for(cn in colnames(cc.matrix.manual)){
         
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
        LOD=do.call('rbind', parLapply(cl, names(thetas), domap) )
        LOD[is.na(LOD)]=0
        LOD=LOD/(2*log(10))
        rownames(LOD)=names(thetas)
        cnn=gsub('/',':', cn)
        saveRDS(LOD, file=paste0(results.dir, 'LOD_NB_', cnn,'.RDS'))
        msig=which(LOD>4, arr.ind=T)
        png(file=paste0(results.dir, 'seg_diagnostic_plot_', cnn, '.png'), width=768,height=768)
        plot(msig[,2], msig[,1])
        dev.off()
    }
}



#joint mapping with additive effect of cell cycle 
for (ee in c(1:4,7:12) ) {
    experiment=experiments[ee]
    print(ee)
    print(experiment)    
    #cross=crossL[[cdesig[ee]]]
    #parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')

    segDataList=readRDS(file=paste0(results.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    #barcode.featues=segDataList$barcode.features

    dispersion.df=segDataList$dispersion.df
    counts=readRDS(paste0(results.dir, 'counts.RDS'))

    # replace if barcode.features is present
    #mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)))
  
    #could switch between thetas here 
    thetas=dispersion.df$theta.cis
    names(thetas)=dispersion.df$gene


    cc.big.table=readr::read_delim('/data/single_cell_eQTL/yeast/results/cell_cycle_v5/cell_cycle_feb02162022.tsv', delim='\t')
    cc.df=cc.big.table %>% dplyr::filter( cell_cycle != "HALPOID" & cell_cycle != "HAPLOIDS" & named_dataset == experiment )
    #cell_cycle != "ALPHA" &
   
    cc.df=cc.df[cc.df$cell_name %in% rownames(Gsub),]
    cc.df$seurat_clusters=as.factor(cc.df$seurat_clusters)
   
    Gsub=Gsub[cc.df$cell_name,]
    Yr=Yr[cc.df$cell_name,]

    print(all.equal(cc.df$cell_name, rownames(Gsub)))

    mmp0=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))



    
    Gsub=segDataList$Gsub[cc.df$cell_name,]
    Yr=segDataList$Yr[cc.df$cell_name,]

    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)
    mmp1=cbind(mmp0, cc.matrix.manual)


    print('calculating LODs')
    #Yks=Matrix(Yr, sparse=T)
    clusterExport(cl, varlist=c("thetas", "Yr", "Gsub", "mmp1", "nbLL", "domap"))
    clusterEvalQ(cl, { Y=Yr;    DM=mmp1;   return(NULL);})
    LOD=do.call('rbind', parLapply(cl, names(thetas), domap) )
    LOD[is.na(LOD)]=0
    LOD=LOD/(2*log(10))
    rownames(LOD)=names(thetas)
    saveRDS(LOD, file=paste0(results.dir, 'LOD_NB_joint.RDS'))
    msig=which(LOD>5, arr.ind=T)
    png(file=paste0(results.dir, 'seg_diagnostic_plot_joint.png'), width=768,height=768)
    plot(msig[,2], msig[,1])
    dev.off()
}





cc.big.table=readr::read_delim('/data/single_cell_eQTL/yeast/results/cell_cycle_v5/cell_cycle_feb02162022.tsv', delim='\t')

#mapping cell-cycle
#for(ee in (7,1:4,8:10,
for (ee in 15:16) {  #1,12,15,16)) {
 #   ee=3
     experiment=experiments[ee]
    print(ee)
    print(experiment)    
    #cross=crossL[[cdesig[ee]]]
    #parents=parentsL[[cdesig[ee]]]
    data.dir=data.dirs[ee]
    results.dir=paste0(base.dir, 'results/', experiment, '/')

    segDataList=readRDS(file=paste0(results.dir, 'segData.RDS'))

    Yr=segDataList$Yr
    Gsub=segDataList$Gsub
    cisMarkers=segDataList$cisMarkers 
    transcript.features=segDataList$transcript.features
    #barcode.featues=segDataList$barcode.features

    dispersion.df=segDataList$dispersion.df
    counts=readRDS(paste0(results.dir, 'counts.RDS'))

    # replace if barcode.features is present
    #mmp1=model.matrix(lm(Yr[,1]~log(barcode.features$nUMI)))
    mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))

    #read cell cycle info
    #cc.df=read_csv(paste0(base.dir, 'results/cell_cycle/', experiment, '/cell_cycle_assignments.csv'))
    cc.df=cc.big.table %>% dplyr::filter(cell_cycle != "ALPHA" & cell_cycle != "HALPOID" & cell_cycle != "HAPLOIDS" & named_dataset == experiment )

    cc.df=cc.df[cc.df$cell_name %in% rownames(Gsub),]
    cc.df$seurat_clusters=as.factor(cc.df$seurat_clusters)
   
    Gsub=Gsub[cc.df$cell_name,]
    Yr=Yr[cc.df$cell_name,]

    print(all.equal(cc.df$cell_name, rownames(Gsub)))

    mmp1=model.matrix(lm(Yr[,1]~log(colSums(counts[,rownames(Yr)]))))
    cc.matrix.manual=with(cc.df, model.matrix(~cell_cycle-1))
    cc.matrix.auto=with(cc.df, model.matrix(~seurat_clusters-1))
    cc.incidence=cbind(cc.matrix.manual, cc.matrix.auto)

    clusterExport(cl, varlist=c("Gsub", "mmp1", "cc.incidence", "domap_logistic"))
    #clusterEvalQ(cl, { Y=Yr;    DM=mmp1;   return(NULL);})
    LOD=do.call('rbind', parLapply(cl, colnames(cc.incidence), domap_logistic) )
    LOD[is.na(LOD)]=0
    LOD=LOD/(2*log(10))
    rownames(LOD)=colnames(cc.incidence)
    dir.create(paste0(base.dir, 'results/cell_cycle_v5/', experiment,'/'))
    saveRDS(LOD, file=paste0(base.dir, 'results/cell_cycle_v5/', experiment, '/cell_cycle_assignment_LOD.RDS'))

    pdf(file=paste0(base.dir, 'results/cell_cycle_v5/', experiment, '/cell_cycle_assignment_LOD.pdf'), width=10, height=5)
    for(i in 1:nrow(LOD)){
        plot(LOD[i,],main=rownames(LOD)[i], ylab='LOD', xlab='marker index')
        abline(v=cumsum(c(0,rle(tstrsplit(colnames(Gsub), '_')[[1]])$lengths)), lty=2, col='blue')
    }
    dev.off()
}





LODo=readRDS(paste0(results.dir, 'LOD_NB3.RDS'))


    #insert code speedup for dispersion estimate here
    L=irlba(Gsub,10)
    colnames(L$u)=paste0('u',1:ncol(L$u))
    mmpN=model.matrix(lm(Yr[,1]~log(colSums(counts))))
    mmpN=cbind(mmpN, L$u)

    pcNB=nebula(counts[expressed.transcripts,], mmpN[,1], pred=mmpN)
    pcModel.df=data.frame(gene=pcNB$summary$gene,theta=1/pcNB$overdispersion[,2])
    nullModel.df=data.frame(gene=pcNBnull$summary$gene,theta=1/pcNBnull$overdispersion[,2])

    cisModel.df=data.frame(gene=names(thetas),theta=thetas)

    test=left_join(cisModel.df, pcModel.df, by='gene', suffix=c('.cis', '.pc'))
    test=left_join(test, nullModel.df, by='gene', suffix=c('','.null'))
    theta_v_expression=left_join(transcript.features, test, by='gene')
    
    library(viridis)
    library(ggpubr)
    a=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta.cis),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()
    b=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta.pc),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()

  
   cc=ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()

    ggarrange(a,b,cc)
plot(log(test$theta.x), log(test$theta.y))

    saveRDS(pcNB, file=paste0(results.dir, 'pcNB3.RDS'))

    pcNBnull=nebula(counts[expressed.transcripts,], mmp1[,1], pred=mmp1)


    test=cv.glmnet(mmpN, Yr[,4], family='poisson', alpha=1)

    
    print('Calculating overdispersions') 
    pcNB=nebula(counts[expressed.transcripts,], mmpN[,1], pred=mmpN)
    saveRDS(pcNB, file=paste0(results.dir, 'pcNB3.RDS'))
    #pcNB=readRDS(paste0(results.dir, 'pcNB2.RDS'))
    thetas=1/pcNB$overdispersion$Cell  
    names(thetas)=pcNB$summary$gene



    theta_v_expression=left_join(pcModel.df, transcript.features, by='gene')
    ggplot(theta_v_expression, aes(x=log2(tot.expression),
                                   y=log2(1/theta),
                                   col=log2(non.zero.cells)))+geom_point()+scale_colour_viridis()
  









       # LODo=readRDS(paste0(results.dir, 'LOD_NB2.RDS'))


}









match(rownames(LOD), sgd.genes$gene_id)


    resids=sapply(cisNB, function(x) x$pearsonResiduals)
    rownames(resids)=rownames(Yr)

    cc.resids=resids[,colnames(resids) %in% cc$ORF]
    um=uwot::umap(cc.resids, n_neighbors=30, metric="cosine", min_dist=0.3, n_threads=36)
    df=data.frame(u1=um[,1], u2=um[,2],data.matrix(Yr))
    ggplot(df,aes(x=u1,y=u2, col=log10(YBL003C)+1))+scale_color_viridis()+geom_point()

    rawC=standardise2(data.matrix(log(Yr[,colnames(Yr) %in% cc$ORF]+1)/log(colSums(counts))))
    um2=uwot::umap(rawC, n_neighbors=30, metric="cosine", min_dist=0.3, n_threads=36)




    df=data.frame(u1=um2[,1], u2=um2[,2],data.matrix(Yr))
    ggplot(df,aes(x=u1,y=u2, col=log10(YBL003C)+1))+scale_color_viridis()+geom_point()





 #   cell.cycle.classification.file=paste0(base.dir, cell.cycle.classification.dir, cell.cycle.classification.names[ee], '.cluster.assignments.txt')
 #   cell.cycle.annot.file=paste0(base.dir, cell.cycle.classification.dir, cell.cycle.classification.names[ee], '.cluster.annotations.txt')
    # get cell cycle covariates
#    cell.covariates=getCCinfo(cell.cycle.classification.file,cell.cycle.annot.file)
    #note n gets smaller here due to cells missing classification
 #   cell.covariates=cell.covariates[!is.na(cell.covariates$cid),]
 #   print(experiment)
 #   print(nrow(cell.covariates))
 #   saveRDS(cell.covariates, paste0(results.dir, 'cell_covariates.RDS'))
  
    
    # par(xaxs='i', yaxs='i')
  #  plot(m.granges$gcoord, colSums(pg[rownames(pg)%in%rownames(vg),])/2727, ylim=c(0,1))
  #  plot(m.granges$gcoord, colSums(vg[classification,]-1)/sum(classification), ylim=c(-1,1), xlab='mpos', ylab='af')
  #  abline(v=gcoord.key, lty=2, col='lightblue')
  #  cis.marker=getCisMarker(sgd.genes, m.granges, counts)
 #vgV=getSavedGenos(chroms, results.dir, type='viterbi')
    #par(xaxs='i', yaxs='i')
    #plot(m.granges$gcoord, colSums(pg[rownames(pg)%in%rownames(vg),])/2727, ylim=c(0,1))
    #plot(m.granges$gcoord, colSums(vgV[classification,]-1)/sum(classification), ylim=c(0,1), xlab='mpos', ylab='af')
    #abline(v=gcoord.key, lty=2, col='lightblue')

     # assuming classification matches experiment, subset and propagate order across phenos and genos 

 #cisNB=fitNBCisModel(cMsubset, mmp1, Yr, Gsub, resids='pearson')          #, resids='deviance')
    #saveRDS(cisNB, file=paste0(results.dir, 'cisNB.RDS'))
    #thetas=sapply(cisNB, function(x) as.numeric(x$theta))
    #intercepts=sapply(cisNB, function(x) as.numeric( x$fmodelBs[1]) )
    #plot(1/comp_models$theta.x, 1/comp_models$theta.y, main='1/theta', xlim=c(0,25), ylim=c(0,25), xlab='cis model', ylab='10 pc model')
    
    #cisModel=readRDS(paste0(results.dir, 'cisNB.RDS'))
    #pcModel=readRDS(paste0(results.dir, 'pcNB2.RDS'))
    #thetas=sapply(cisModel, function(x) as.numeric(x$theta))
    #cisModel.df=data.frame(gene=names(thetas), theta=thetas)
    #pcModel.df=data.frame(gene=pcModel$summary$gene,theta=1/pcModel$overdispersion[,2])
    #comp_models=left_join(cisModel.df, pcModel.df, by='gene')











    ## optional get ASE info
    # ase.counts=getASEcounts(g.counts, sgd.genes, m.granges, 16)
    # tot.geno.inf=ase.counts[[1]]+ase.counts[[2]] 

    # taf=colSums(ase.counts$ref)/(colSums(ase.counts$ref+ase.counts$alt))
    # caf=log2(colSums(ase.counts$ref+ase.counts$alt))
    # plot(sgd.genes$gcoord[match(names(taf), sgd.genes$ID)],taf, col=(caf<5)+1)
    # abline(v=gcoord.key, lty=2, col='lightblue')
    # hist(log2(colMeans(tot.geno.inf)), xlab='log2(avg genotype informative counts per transcript per cell)')


    #Gcor=crossprod(G)/(nrow(G)-1)
    #rm(G)
    # visualize auto-correlation. e.g. correlation decay from given marker to a marker with an index 100 less
    #v=Gcor[row(Gcor)==(col(Gcor)-100)] ... eek
    # prune correlated markers  (approx every 5cm)
    # previous code transformed genotypes as log(g/(1-g)) for a beta-regression-like setup
    # ~30 or markers were lost by div by 0 errors, fixed 2/14
    #fcm=caret::findCorrelation(Gcor,cutoff=.999, verbose=F)














 
    ## get genotyypes
    #vg=getSavedGenos(chroms, results.dir, type='viterbi')
    #saveRDS(vg, paste0(results.dir, 'vit_gmatrix.RDS'))

    #pg=getSavedGenos(chroms, results.dir, type='genoprobs')
    
        # antiquated, read in r/qtl binned formatted genotypes   
    #rQTL.coded.file=paste0(results.dir, 'rQTLcoded.RDS')
    #rQTL.coded=readRDS(rQTL.coded.file)


















    #plot(m.granges$gcoord, colSums(vg-1)/nrow(vg), ylim=c(0,1))
    #abline(v=gcoord.key)

    cell.covariates$total.counts=colSums(counts)
    
    Yr=t(as.matrix(counts))
    Yr.var=colVars(Yr, parallel=T)
    Yr=Yr[,Yr.var>0]
    tcounts=Yr>0
    tcounts=colSums(tcounts)
    Yr=Yr[,tcounts>20]

    countsR=counts[colnames(Yr),]

    #.05 230
    mm=model.matrix(counts[1,]~log2(total.counts)+cid,data=cell.covariates) 
    #Yre2=matrix(NA, nrow(Yr), ncol(Yr))
    #colnames(Yre2)=colnames(Yr)
    #rownames(Yre2)=rownames(Yr)
    #for(g in 1:ncol(Yr)){
    #    print(g)
    #    Yre2[,g]=scale(residuals(glm(Yr[,g]~offset(log(total.counts))+cid,data=cell.covariates, family=poisson()),'pearson' ))
    #}
    
    #mm0=model.matrix(counts[1,]~log2(total.counts),data=cell.covariates) 
    #BOLS=lm.fit(mm0, log2(Yr+1))
    #Yr0=scale(residuals(BOLS))
    #rPC=hd.eigen(t(Yr0),vectors=T, center=F, scale=F, k=2)
    #plot(rPC$vectors[,1], rPC$vectors[,2])

    BOLS=lm.fit(mm, log2(Yr+1))
    Yresid=residuals(BOLS)
    rownames(Yresid)=rownames(Yr)
    colnames(Yresid)=colnames(Yr)
    rm(BOLS)

    Yre=standardise2(Yresid)
    rownames(Yre)=rownames(Yresid)
    colnames(Yre)=colnames(Yresid)
    #------------------------------------------------------------------------------------------------------------

    #Yre=Yre2
    assignme=model.matrix(counts[1,]~cell.covariates$cid-1)
    cgen=scale(residuals(lm(assignme~log2(cell.covariates$total.counts))))

    rcc=crossprod(cgen,G)/(nrow(Yre)-1)
    LODcc=-nrow(Yre)*log(1-rcc^2)/(2*log(10))
    #saveRDS(LODcc, file=paste0(results.dir, 'lodCC.RDS'))
    par(mfrow=c(7,1), xaxs='i')
    for(i in 1:7) {
    plot(m.granges$gcoord, LODcc[i,], main=rownames(rcc)[i], ylab='LOD')
    abline(v=gcoord.key, lty=2, col='lightblue')
    }


   # 1D scan speedup 
    r=crossprod(Yre,G)/(nrow(Yre)-1)
    LOD=-nrow(Yre)*log(1-r^2)/(2*log(10))

    #pool for permutations
    wgFDR= getFDRfx(r, Yre, G, nperm=5)
    peaks1D=get1Dpeaks(chroms,r,LOD,wgFDR, m.granges, sgd.genes)
    saveRDS(peaks1D, file=paste0(results.dir, 'peaks1D.RDS'))
    plot2D(peaks1D, .1, gcoord.key=gcoord.key, experiment=experiment)
   
   
    #multivariate scan for hotspots 
    mLOD=list()
    for(cc in chroms) {
        print(cc)
        ng=sgd.genes[seqnames(sgd.genes)!=cc,]
        YreS=Yre[,(colnames(Yre) %in% ng$Name)]
        moi=m.granges$sname[as.character(seqnames(m.granges))==cc]
        #test=hd.eigen(t(Yre), center=F, scale=F,vectors=T)
        pc.to.retain=40
        yreduced=hd.eigen(t(YreS), center=F, scale=F, k=pc.to.retain, vectors=T)
        testF=mvn.scanone(G[,moi], yreduced$vectors)
        testN=determinant(crossprod(yreduced$vectors), logarithm=T)$modulus
        mLODv=(nrow(yreduced$vectors)/2)*(testN-testF)/(2*log(10))
        mLOD[[cc]]=mLODv
    }
    mLOD=do.call('c', mLOD)    
    peaksMV=getMVpeaks(chroms, r, LOD, mLOD, wgFDR, m.granges, sgd.genes)
    #saveRDS(peaksMV, file=paste0(results.dir, 'peaksMV.RDS'))


    ccPeaks=list()
    # split by cell -cycle
    cell.covs=split(cell.covariates, cell.covariates$cid)
    #Yr and vg are counts and genos respectively
    for(covn in names(cell.covs)){
            print(covn)
          ccdf=cell.covs[[covn]]
          Yrs=Yr[ccdf$barcode,]
          tcounts=Yrs>0
          tcounts=colSums(tcounts)
          Yrs=Yrs[,tcounts>20]
          
          mms=model.matrix(Yrs[,1]~log2(total.counts),data=ccdf) 
          BOLSs=lm.fit(mms, log2(Yrs+1))
          Yresids=residuals(BOLSs)
          rm(BOLSs)
          rownames(Yresids)=rownames(Yrs)
          colnames(Yresids)=colnames(Yrs)
          Yres=standardise2(Yresids)
          rownames(Yres)=rownames(Yresids)
          colnames(Yres)=colnames(Yresids)
        
          vgs=vg[ccdf$barcode,]
          Gs=standardise2(vgs)
          rs=crossprod(Yres,Gs)/(nrow(Yres)-1)
          LODs=-nrow(Yres)*log(1-rs^2)/(2*log(10))

         #pool for permutations
         wgFDRs= getFDRfx(rs, Yres, Gs, nperm=5)
         peaks1Ds=get1Dpeaks(chroms,rs,LODs,wgFDRs, m.granges, sgd.genes)
         ccPeaks[[covn]]=peaks1Ds
    }
    saveRDS(ccPeaks, file=paste0(results.dir, 'CCpeaks.RDS'))



}
texp=colSums(Yr)
llr=apply(ll,1,range)
llm=apply(ll,1,mean)
plot(log10(texp), max.obsLOD,col='red', ylim=c(0, 0.2),xlim=c(1,6.5),
     xlab='log10(total expression per transcript)',
     ylab='abs(r)'
    )
points(log10(texp), llm, ylim=c(0.02,0.06))
segments(log10(texp),llr[1,],log10(texp),llr[2,], col="#00000066")
points(log10(texp), max.obsLOD,col='red', ylim=c(0, 0.2))
#points(log10(texp), llm,  col="#00000022")

qm=rbindlist(ccPeaks, idcol='cell_cycle')
qm=qm[qm$FDR<.1,]
ggplot(qm,aes(x=mgcoord,y=tgcoord, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    geom_segment(aes(x =mgcoordL, y = tgcoord, xend =mgcoordR, yend = tgcoord)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    facet_wrap(~cell_cycle)+
    geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/100919/yeast_CC_QTL.png')


countsR=as.matrix(countsR)

p1dr=peaks1D[peaks1D$FDR<.1,]
p1dr$transcript=droplevels(p1dr$transcript)
gqtl=split(p1dr, p1dr$transcript)

cell.covariates$cid=relevel(cell.covariates$cid, "G1")

cc.fits=foreach(gene=names(gqtl)) %do% {
    # cell cycle 
    tryCatch({
    #gene=names(gqtl)[3]

    y=countsR[gene,]
    XX=G[,as.character(gqtl[[gene]]$peak.marker)]
    #XX=matrix(XX)
    #ZXX=Z%*%XX
    #colnames(XX)=as.character(gqtl[[gene]]$peak.marker) #sb[[gene]]$fscan.marker
    yp=log2(y+1)

    obs=data.frame(y=y, yp=yp, cell_total=cell.covariates$total.counts, Cid=cell.covariates$cid)

    #m2=model.matrix(y~cell.covariates$cid-1)
    #mm0=lm(yp~log2(cell_total)+XX,data=obs)
    #mm0=lm(yp~log2(cell_total)+Cid+XX,data=obs)

    #given this is multithreaded, do not use dopar!
    mm=lm(yp~log2(cell_total)+XX*Cid,data=obs)
    #mm=glm.nb(y~offset(log(cell_total))+XX*XX*Cid,data=obs)
    #mm=glm.nb(y~(log(cell_total))+XX*Cid,data=obs)

    clist=tidy(mm)
    if(is.matrix(XX)==FALSE) {
      clist$term=gsub('XX',   paste0('XX', as.character(gqtl[[gene]]$peak.marker)), clist$term)
    }
    print(gene)
    print(data.frame(clist))
     return(clist)
    },error=function(e) {
        return(NULL)
    })
}
names(cc.fits)=names(gqtl)
saveRDS(cc.fits, file = paste0(results.dir, 'CCfits.RDS'))
#save(cc.fits, file='/data/single_cell_eQTL/yeast/cc.fits.RData')
# for regular lm
cc.ps=sapply(cc.fits, function(x) {
                 y=x$p.value[grep(':', x$term)] 
                 names(y)=x$term[grep(':', x$term)]
                 return(y)
    } )
hist(unlist(cc.ps),breaks=350, xlab='p', main='(QTL x Cell_Cycle) p-values')

rcc.ps=lapply(cc.ps,function(n) {
           ss=strsplit(names(n),':')
           qm=sapply(ss, function(x) gsub('XX', '', x[1]))
           cm=sapply(ss, function(x) x[2])
           data.frame(peak.marker=qm,cell_cycle=cm,p=as.vector(n),stringsAsFactors=F)
    })
rcc.ps=rbindlist(rcc.ps, idcol='transcript')
rcc.ps$FDR=qvalue(rcc.ps$p)$qvalue

    rcc.ps$mchr=as.character(seqnames(m.granges))[match(rcc.ps$peak.marker, m.granges$sname)]
    rcc.ps$mpos=(start(m.granges))[match(rcc.ps$peak.marker, m.granges$sname)]
    rcc.ps$mgcoord =m.granges$gcoord[match(rcc.ps$peak.marker, m.granges$sname)]
    #rcc.ps$mgcoordL=m.granges$gcoord[match(rcc.ps$CI.l, m.granges$sname)]
    #rcc.ps$mgcoordR=m.granges$gcoord[match(rcc.ps$CI.r, m.granges$sname)]
    #rcc.ps$mposL   =(start(m.granges))[match(rcc.ps$CI.l, m.granges$sname)]
    #rcc.ps$mposR   =(start(m.granges))[match(rcc.ps$CI.r, m.granges$sname)]

    rcc.ps$tchr=as.character(seqnames(sgd.genes))[match(rcc.ps$transcript, sgd.genes$Name)]
    rcc.ps$tpos=start(sgd.genes)[match(rcc.ps$transcript, sgd.genes$Name)]
    rcc.ps$tgcoord=sgd.genes$gcoord[match(rcc.ps$transcript, sgd.genes$Name)]
saveRDS(rcc.ps, file=paste0(results.dir, 'CCpeaksxQTL_int.RDS'))

qm2=rcc.ps[rcc.ps$FDR<.1,]

ggplot(qm2,aes(x=mgcoord,y=tgcoord, alpha=(-log10(FDR+1e-6)/6)))+geom_point()+
    #geom_segment(aes(x =mgcoordL, y = tgcoord, xend =mgcoordR, yend = tgcoord)) + #, color ="grey70", alpha=.1))+
    xlab('')+ylab('')+ scale_alpha(guide = 'none')+
    facet_wrap(~cell_cycle)+
    geom_hline(yintercept=gcoord.key,color='lightblue')+
    geom_vline(xintercept=gcoord.key,color='lightblue')+
    theme_classic()+
   theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())           

sq=split(qm2, qm2$transcript)
sq=sq[which(sapply(sq,nrow)>1)]
sqf=lapply(sq,function(x) x[x$mchr!='chrIV',])
t1=cc.fits[names(which(sapply(sqf,nrow)>1))]

#svF=irlba(G, nv=150)
#svFGG=irlba(G*G, nv=150)

disp.fits=foreach(gene=names(gqtl)) %dopar% {

 tryCatch({
    #gene=rownames(counts)[r]
    #  dispersion test 
    #gene='YHB1'
    #gene='YJR009C'
    print(gene)
    #y=raw_data[r,]
    y=countsR[gene,]
    XX=G[,as.character(gqtl[[gene]]$peak.marker)]
    yp=log2(y+1)
    obs=data.frame(y=y, yp=yp, cell_total=cell.covariates$total.counts, Cid=cell.covariates$cid)
    #test=model.matrix(y~(XX[,1]+XX[,2]+XX[,3]+XX[,4])^2-1)
    bm=glmmTMB(y~offset(log(cell_total))+Cid+XX, family=nbinom2,data=obs)
    #bm2=glmmTMB(y~offset(log(cell_total))+Cid+test, family=nbinom2,data=obs)

    # bm=glmmTMB(y~(log(cell_total))+Cid+XX, family=nbinom2,data=obs)
    #bm2=update(bm, dispformula=~Cid+XX)


   if(is.null(dim(XX)[2])) {XX=as.matrix(XX)}
   dlist=list()
   for(i in 1:ncol(XX)){
            bm1=update(bm, dispformula=~XX[,i])
            ps=as.numeric(pchisq(2*(logLik(bm1)-logLik(bm)),1,lower.tail=F))
            #nvec=c(summary(bm1)$coefficients$cond[i+4,], sigma(bm1), ps)
            nvec=c(summary(bm1)$coefficients$cond[i+4,], sigma(bm1),ps,summary(bm1)$coefficients$disp[2,])
            names(nvec)[5]='theta_int'
            names(nvec)[6]='theta_coef'
            names(nvec)[7]='lrt_p'
            names(nvec)[8:11]=paste0('theta_coef_',names(nvec)[8:11]) 
            dlist[[as.character(gqtl[[gene]]$peak.marker)[i]]]=nvec
        }

    print(dlist)
    return(dlist)

    },error=function(e) {
        dlist=list()
        nvec=rep(NA,11)
        names(nvec)=c("Estimate",  "Std. Error",            "z value",               "Pr(>|z|)"  ,           
         "theta_int"  ,           "theta_coef"      ,      "lrt_p"   ,              "theta_coef_Estimate" , 
         "theta_coef_Std. Error", "theta_coef_z value"  ,  "theta_coef_Pr(>|z|)")  
        dlist[['fail']]=nvec
        return(dlist)
    })
}
names(disp.fits)=names(gqtl)
hist(unlist(sapply(disp.fits, function(x) sapply(x, function(y)y[4]))),breaks=100)
hist(unlist(sapply(disp.fits, function(x) sapply(x, function(y)y[7]))),breaks=100)


dp3df=rbindlist(lapply(disp.fits, function(x) data.frame(t(x[[1]]))),idcol='gene')

df.table=rbindlist(disp.fits,idcol='transcript',fill=T)

df.table=rbindlist(disp.fits,idcol='transcript')




    bm=glmmTMB(y~offset(log(cell_total))+Cid*XX, family=nbinom2,data=obs)
    clist=tidy(bm)
    bm2=glm.nb(y~offset(log(cell_total))+Cid*XX,data=obs)
    
    #counts[1,]~log2(total.counts)+cid,data=cell.covariates) 
    



h=t2$cPeaks[t2$cPeaks$chr=='chrVIII',]
par(xaxs='i')
plot(m.granges$gcoord, mLOD, type='l')
abline(v=gcoord.key, lty=2, col='lightblue')

     
# chrVIII is gpa1
# ylim=c(0,30), type='l')  # /4.6))


plot2D=function(peaks1D, FDR.thresh=.05, gcoord.key=gcoord.key, experiment=experiment) {
    plot.thresh=FDR.thresh
    par(xaxs='i', yaxs='i')
    plot(peaks1D$mgcoord[peaks1D$FDR<plot.thresh], peaks1D$tgcoord[peaks1D$FDR<plot.thresh],
         ylab='transcript position',
         xlab='QTL position', main=experiment, xaxt='n', yaxt='n', pch=20)
    axis(1, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    axis(2, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    abline(v=gcoord.key, lty=2, col='grey')
    abline(h=gcoord.key, lty=2, col='grey')
}







peaks1D=readRDS('/data/single_cell_eQTL/yeast/results/2_2444_44-/peaks1D.RDS')


 


















    # a chromosome of interest for plotting 
    coii='chrVII'

    pp.file=paste0(results.dir, coii,'.RDS')
    post.prob=readRDS(pp.file)
    
    v.file=paste0(results.dir, coii, '_viterbi.RDS')
    vit= readRDS(v.file)

    for(i in 100:200) {
        diagnoseHMM(i,chrom=coii, gmap.ss, g.counts[[1]], g.counts[[2]], post.prob, rQTL.coded, viterbiPath=vit, classification=classification)
        readline()
    }

}
















    # recurring.ref=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==1,na.rm=T))
    # recurring.alt=apply(rQTL.coded[,-het.cells],1,function(x) sum(x==2,na.rm=T))
    # raf=recurring.ref/(recurring.alt+recurring.ref)
    # taf=recurring.ref+recurring.alt
    # approx rate of hets vs homzyg = 0.0075
    # 29234/(2016111 +1881594)



    #saveRDS(classifications, file = paste0(base.dir, 'is_haploid_classifcation.RDS'))
#################################################

 
# usually we'll want to parallelize this further

#v=viterbi(cross2, error_prob=.01, ,lowmem=F,quiet=F, cores=16)
#dv=do.call('cbind', v)
#by.freq=apply(dv,2, function(x) sum(x==1)/length(x))
#dv.prune=dv[,by.freq<.7 & by.freq>.3]
#G=dv.prune

G2=do.call('cbind', sapply(gps, function(x) x[,2,]))

YJMpS=rowSums(G2)


Ys=t(as.matrix(counts))
Ys=Ys[-het.cells,]

total.counts=rowSums(Ys)

mm=model.matrix(Ys[,1]~log2(total.counts))
BOLS=lm.fit(mm, log2(Ys+1), intercept=F)
Yresid=residuals(BOLS)
rownames(Yresid)=rownames(Ys)
colnames(Yresid)=colnames(Ys)

Yre=standardise(Yresid)
#Yre=standardise(Rfast::colRanks(Yresid))
rownames(Yre)=rownames(Yresid)
colnames(Yre)=colnames(Yresid)
nac=apply(Yre, 2, function(x) sum(is.na(x)))
Yre=Yre[,nac==0]


Gr=standardise(G2)
colnames(Gr)=colnames(G2)
rownames(Gr)=rownames(Yre)

r=crossprod(Yre,Gr)/(nrow(Yre)-1)
LOD=(-nrow(Yre)*log(1-r^2))/(2*log(10))
sigL=which(LOD>8, arr.ind=T)
plot(sigL[,2], sigL[,1], pch=21)

getFDRfx=function(r, Yre,Gr,nperm=5, vint=seq(0.001,.15,.001)){
    
    max.obsLOD=apply(abs(r),1,max)

    permR=replicate(nperm, {
                    nind=sample(1:nrow(Yre))
                    crossprod(Yre[nind,],Gr)/(nrow(Yre)-1)
    })

    ll=apply(abs(permR),3, function(x) apply(x,1,max))
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = vint

    # expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(apply(ll,2,function(x) sum(x>thresh))))
                 })
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    

    pFDR = rev(cummax(rev(pFDR)))
    fdrFX=approxfun(pFDR, vint, ties=mean)
    rtoFDRfx=approxfun(vint,pFDR, ties=mean)

    return(list(fdrFX=fdrFX,rtoFDRfx=rtoFDRfx))
}
#getPeaks=function(r, thresh){
#    cind=gsub('_.*','', colnames(r))
#    cmaxind=matrix(0, nrow(r), 16)
#    rownames(cmaxind)=rownames(r)
#    cmax=cmaxind
#    for(i in 1:nrow(r)) {
#        #print(i)
#        x=r[i,]
#      y=split(x, cind);
#     cmaxind[i,]=as.vector(sapply(y, function(h) names(which.max(abs(h)))))
#     cmax[i,]=x[cmaxind[i,]]
#    }
#    tmax=rep(rownames(cmaxind),16)
#    sig.hits=which(abs(cmax)>thresh) #,arr.ind=T)
#
#    qtl.p=cmaxind[sig.hits]
#    qtl.t=tmax[sig.hits]
#    return(data.frame(transcript=qtl.t,marker=qtl.p, stat=cmax[sig.hits]))
#}

#Lthresh=(nrow(Yre)*log(1-ff(.05))^2)/(2*log(10))
chromFDR=sapply(chroms, function(x) {
                 print(x)
                 Gs=Gr[,grep(paste0('^',x,'_'), colnames(Gr))]
                 r=crossprod(Yre,Gs)/(nrow(Yre)-1)
                 getFDRfx(r, Yre, Gs, nperm=5)
                 })


cPeaks=list()
for(cc in chroms) {
    print(cc)
    rS=r[,grep(paste0('^',cc,'_'), colnames(r))]
    mstat=apply(rS,1,max)
    mstat.pos=apply(rS,1,which.max)
    cPeaks[[cc]]=data.frame(
    transcript=rownames(rS),
    marker=colnames(rS)[mstat.pos],
    stat=mstat,
    FDR=chromFDR[,cc][[2]](mstat)
    )
    cPeaks[[cc]]$FDR[is.na(cPeaks[[cc]]$FDR)]=1
}


cP=do.call('rbind', cPeaks)
plot(match(cP$marker[cP$FDR<.1], colnames(Gr)),
     match(cP$transcript[cP$FDR<.1], colnames(Yre) ), 
             xlab='marker index', ylab='transcript index', main='joint analysis FDR < 10% corrected counts')

#viterbi = 894 at fdr<10%
#experiments=c('1_2444_44_1-2', '2_2444_44-', '3_ByxRM_51-_1-2', 
#              '4_ByxRM_51-', '0_BYxRM_480MatA_1', '0_BYxRM_480MatA_2')
#cdesig=c('B','B', 'A', 'A', 'A', 'A')
#cell.cycle.classification.names=c('set1', 'set2', 'set3', 'set4', 'set01', 'set02')
#cell.cycle.classification.dir='results/cell_cycle/yeast_cc_annotations_v5/'


#experiments=c('08_2444_cross_10k_Feb_21',
#              '09_2444_cross_5k_Feb_21',
#              '10_3004_cross_10k_Feb_21',
#              '11_3004_cross_5k_Feb_21',
#              '07_2444-cross-1',
#              '07_2444-cross-2'
#        )
#cdesig=c('B', 'B',  '3004', '3004', 'B', 'B')           
#              '1_2444_44_1-2', '2_2444_44-', '3_ByxRM_51-_1-2', 
#              '4_ByxRM_51-', '0_BYxRM_480MatA_1', '0_BYxRM_480MatA_2')
#cdesig=c('B','B', 'A', 'A', 'A', 'A')
#cell.cycle.classification.names=c('set1', 'set2', 'set3', 'set4', 'set01', 'set02')
#cell.cycle.classification.dir='results/cell_cycle/yeast_cc_annotations_v5/'
#by visual inspection (might want to revist this classification)
                # for 08-11,07
#het.thresholds=c(2^5.2, 2^4.9, 2^5.6, 2^4.9, 2^6,2^5.5)
                #for 1,2,3,4,0,0
                 #c(2^6.1, 2^7.5, 2^7, 2^7.4, 2^6.5,2^6.5)


#-----------------------------------------------------------------------------


