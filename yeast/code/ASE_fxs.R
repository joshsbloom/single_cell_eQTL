library(GenomicRanges)
library(rtracklayer)
library(seqinr)
library(spam)
library(data.table)
library(vcfR)

crosses.to.parents=list(
     '375'=c("M22", "BYa"),
     'A'  =c("BYa", "RMx"),
     '376'=c("RMx", "YPS163a"),
     'B'  =c("YPS163a", "YJM145x"),
     '377'=c("YJM145x", "CLIB413a"),
     '393'=c("CLIB413a", "YJM978x"),
     '381'=c("YJM978x", "YJM454a"),
    '3008'=c("YJM454a", "YPS1009x"),
    '2999'=c("YPS1009x", "I14a"),
    '3000'=c("I14a", "Y10x"),
    '3001'=c("Y10x", "PW5a"),
    '3049'=c("PW5a", "273614xa"),
    '3003'=c("273614xa", "YJM981x"),
    '3004'=c("YJM981x", "CBS2888a"),
    '3043'=c("CBS2888a", "CLIB219x"),
    '3028'=c("CLIB219x", "M22")
    )
chroms=paste0('chr', as.roman(1:16)) 

reference.dir='/data/single_cell_eQTL/yeast/reference/'
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

#c('I','II', 'III', 'IV', 'V', 'X')

# Get Gene Intervals

#Previous code 
#sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
#reference.dir='/data/single_cell_eQTL/yeast/reference/'
#sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
#sgd=as.data.frame(sgd.granges)
#gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

#sgd.genes=sgd.granges[sgd.granges$type=='gene',]
#sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)

#update with new gff to contain 3' UTRs
#sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
#sgd=as.data.frame(sgd.granges)
#gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))

#sgd.genes=sgd.granges[sgd.granges$type=='gene',]
#sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)
#--------------------------------------------------------------------------------
getSGD_GeneIntervals=function(reference.dir){
    #update gene definitions with utr
    sgd.granges=import.gff(paste0(reference.dir, 'genes.gtf'))
    sgd=as.data.frame(sgd.granges)
    sgd.genes=sgd.granges[sgd.granges$type=='exon',]
    sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)
    return(sgd.genes)
}



getMarkerIntervals=function(variants) {
    markerGR=makeGRangesFromDataFrame(data.frame(chr=variants[[1]], 
                                                 start=as.numeric(variants[[2]]),
                                                 end=as.numeric(variants[[2]]), 
                                                 strand="*"))
    markerGR$name=paste0(variants[[1]], ':', variants[[2]]) 
    markerGR$gcoord=gcoord.key[as.character(seqnames(markerGR))]+start(markerGR)
    markerGR$name2=paste0(variants[[1]], '_', variants[[2]])
    return(markerGR)
}

buildASE_data=function(data.dir,experiment.name){
    barcodes=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/barcodes.tsv'),sep='\t',header=F,stringsAsFactors=F)[,1]
    features=read.csv(paste0(data.dir,'filtered_feature_bc_matrix/features.tsv.gz'), sep='\t',header=F, stringsAsFactors=F)

    counts=spam::read.MM(paste0(data.dir,'filtered_feature_bc_matrix/matrix.mtx.gz')) #,sep='\t',header=F,stringsAsFactors=F)[,1]
    #umis per cell
    numi=colSums(counts)


    #test that we can run cleanup at this step
    aC=cleanup(spam::read.MM(paste0(data.dir, 'alt_counts.mtx')))
    rC=cleanup(spam::read.MM(paste0(data.dir, 'ref_counts.mtx')))
    #aC=readMM(paste0(data.dir, 'alt_counts.mtx'))
    #rC=readMM(paste0(data.dir, 'ref_counts.mtx'))

    # one idea is to add a filter here to get rid of cells with low overall depth 
    variants=read.csv(paste0(data.dir,'out_var.txt'),sep='\t',header=F,stringsAsFactors=F)[,1]
    variants= tstrsplit(variants,'_', type.convert=T)
    variants[[2]]=variants[[2]]+1

    markerGR=getMarkerIntervals(variants)

    # modified this to only look for expected variants
    #all.vcf=read.vcfR(paste0(reference.dir, 'parents.nostar.vcf'))
    all.vcf=read.vcfR(paste0(data.dir, 'parents.vcf'))
    gts=extract.gt(all.vcf)
    rm(all.vcf)

    #par.vcf=read.vcfR(paste0(reference.dir, 'BYxRMxYPS163xYJM145.vcf'))
    #gts=extract.gt(par.vcf)

    #gts=gtsp[rownames(gts),]#c(2,9,11,16)]
    gts[is.na(gts)]="0"
    gts.subset=gts[rownames(gts)%in%markerGR$name2,]
    
    bad.markers=which(!markerGR$name2%in%rownames(gts.subset))
    if(length(bad.markers)>0) {
    #careful, assuming bad.markers is not empty
        markerGR=markerGR[-bad.markers,]
        aC=aC[-bad.markers,]
        rC=rC[-bad.markers,]
        gts=gts.subset
    }

    um=read.csv(paste0(data.dir, 'analysis/umap/2_components/projection.csv'))
    names(um)[1]='barcode'
    return(list(barcodes=barcodes,
                features=features,
                counts=counts,
                numi=numi,
                um=um,
                aC=aC,
                rC=rC,
                variants=variants,
                markerGR=markerGR,
                gts=gts,
                experiment.factor=rep(experiment.name,length(barcodes))))
}


getLL=function(rC, aC, p1.ref, vname, ee=.002){
    #ivec=gts[,1]=='0'
    #p1.ref=gts[,1]=='0'

    p1=rbind(rC[p1.ref,], aC[!p1.ref,])
    p1=as.dgCMatrix.spam(p1)
    vscramb=c(vname[p1.ref], vname[!p1.ref])
    rownames(p1)=vscramb
    p1=p1[vname,]
    p1=cleanup(as.spam.dgCMatrix(p1))

    p2=rbind(rC[!p1.ref,],aC[p1.ref,])
    p2=as.dgCMatrix.spam(p2)
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


parents.to.crosses=data.frame(
                              cross=rep(names(crosses.to.parents), sapply(crosses.to.parents, length)),
                              parent=unlist(crosses.to.parents), stringsAsFactors=F, row.names=NULL)

#input.dips=table(parents.to.crosses[parents.to.crosses$parent %in% colnames(ase.Data$gts),]$cross)
#input.dips=names(input.dips)[input.dips>1]


getDipAssignments=function(ase.Data, input.diploids, crosses.to.parents, ncores=8){
    registerDoMC(cores=ncores)

    lls=foreach(nn =colnames(ase.Data$gts) ) %dopar% {
         print(nn)
         getLL(ase.Data$rC,ase.Data$aC, ase.Data$gts[,nn]=='0',ase.Data$markerGR$name2)
    }

    do.call('cbind', lls)->llik.table
    colnames(llik.table)=colnames(ase.Data$gts)

    dipLik=list()
    for(dip in input.diploids)  {
        print(dip)
        dipLik[[dip]]=rowSums(llik.table[,crosses.to.parents[[dip]]])
    }

    dipLik=do.call('cbind', dipLik)
   
    m2diff=apply(dipLik,1, function(x) {
        y=sort(x,decreasing=T)
        (y[1]-y[2])
    })

    mcountl=apply(dipLik,1, which.max)
    
    colnames(dipLik)=paste0('llik.',colnames(dipLik))

    df=data.frame(barcode=ase.Data$barcode, 
                  umis_per_cell=ase.Data$numi,
                  dipLik,
                  diploid_assignment_likdiff=m2diff,
                  diploid_assignment=mcountl,
                  diploid = as.vector(sapply(crosses.to.parents[gsub('llik.','',colnames(dipLik))], paste, collapse=' x '))[mcountl],
                  diploid_name= gsub('llik.','',colnames(dipLik))[mcountl]
                  )
    return(df)
}




#some assumptions here are that parents 
countDiploidSpecificVariants=function(ase.Data,input.diploids,crosses.to.parents){
    tcount=colSums(ase.Data$count)
    acount=apply(ase.Data$gts, 1, function(x) sum(x=='1'))
    rvars=acount==1
    gts.sub=ase.Data$gts[rvars,]
    ap=list()
    for(dip in input.diploids) { 
        ap[[dip]]=(gts.sub[,crosses.to.parents[[dip]][1]]=="1" | 
        gts.sub[,crosses.to.parents[[dip]][2]]=="1")
    }
    #ap2=(gts.sub[,3]=="1" | gts.sub[,4]=="1")
    #ap3=(gts.sub[,5]=="1" | gts.sub[,6]=="1")

    aCr=ase.Data$aC[rvars,]
    
    apc=lapply(ap,function(x) colSums(aCr[which(x),]>0))

    df=data.frame(barcode=ase.Data$barcode, do.call('cbind',apc), stringsAsFactors=F)
    names(df)[2:ncol(df)]=paste0('rcount.',input.diploids)
    return(df)
}



getPhasedCountsPerTranscript=function(ase.Data, dip.Assignments, dip,
                                      sgd.genes, crosses.to.parents, threads=32) {

    classified.cells=dip.Assignments$diploid_name==dip
    selected.barcodes=dip.Assignments$barcode[dip.Assignments$diploid_name==dip]
    #phase ASE data
    rCs=as.dgCMatrix.spam(ase.Data$rC)[,classified.cells]
    aCs=as.dgCMatrix.spam(ase.Data$aC)[,classified.cells]

    #gts.subset=gts[rownames(gts)%in%markerGR$name2,]

    # rejigger logic , 5/10/21
    #cross.index=2
    #cross.index=14
    #crosses.to.parents[[cross.index]]
    p1.ref=ase.Data$gts[, crosses.to.parents[[dip]][1]]=='0' 
    p1.not.p2=ase.Data$gts[, crosses.to.parents[[dip]][2]]!=ase.Data$gts[, crosses.to.parents[[dip]][1]]

    p1.index.ref=p1.ref & p1.not.p2
    p1.index.nonref=!p1.ref & p1.not.p2 

    vname=ase.Data$markerGR$name2

    p1=rbind(rCs[p1.index.ref,], aCs[p1.index.nonref,])
    #p1=as.dgCMatrix.spam(p1)
    vscramb=c(vname[p1.index.ref], vname[p1.index.nonref])
    rownames(p1)=vscramb
    #reorder
    p1=p1[vname[(vname %in% vscramb)],]
    #p1=cleanup(as.spam.dgCMatrix(p1))

    p2=rbind(rCs[p1.index.nonref,],aCs[p1.index.ref,])
    #p2=as.dgCMatrix.spam(p2)
    vscramb=c(vname[p1.index.nonref], vname[p1.index.ref])
    rownames(p2)=vscramb
    p2=p2[vname[(vname %in% vscramb)],]
    #p2=cleanup(as.spam.dgCMatrix(p2))

    g.counts=list(p1,p2)

    m.granges=ase.Data$markerGR[ase.Data$markerGR$name2 %in% rownames(p1)]


    fO=findOverlaps(sgd.genes, m.granges)
    #mgsplit=split(m.granges$name[subjectHits(fO)], sgd.genes$Name[queryHits(fO)])
    mgsplit=split(m.granges$name2[subjectHits(fO)], sgd.genes$gene_id[queryHits(fO)])

    print('calc counts per gene, allele 1')

    # could restructure these as sparse matrices if necessary
    par1.ASEcounts=mcmapply( function(x) {
                #m=match(x, markerGR$name2)                               
                if(length(x)>1) {
                    colSums(g.counts[[1]][x,])
                } else {
                    g.counts[[1]][x,]
                }
           },mgsplit,mc.cores=threads)
    rownames(par1.ASEcounts)=selected.barcodes
    print('calc counts per gene, allele 2')
    par2.ASEcounts=mcmapply( function(x) {
                #m=match(x, markerGR$name2)                               
                if(length(x)>1) {
                    colSums(g.counts[[2]][x,])
                } else {
                    g.counts[[2]][x,]
                }
    }, mgsplit,mc.cores=threads)
    rownames(par2.ASEcounts)=selected.barcodes

    return(list(par1.ASEcounts=par1.ASEcounts, par2.ASEcounts=par2.ASEcounts))
}





doBetaBinomialTest=function(dip,phasedCounts,dip.Assignments,ase.Data,
                            nUMI.thresh=10000,informativeCellThresh=64,threads=48) {
  #nUMI.thresh=10000
  #informativeCellThresh=64
  classified.cells=dip.Assignments$diploid_name==dip
  selected.barcodes=dip.Assignments$barcode[dip.Assignments$diploid_name==dip]
  cells.to.keep=ase.Data$numi[classified.cells]<nUMI.thresh 
  p1=phasedCounts[[1]][cells.to.keep,]
  p2=phasedCounts[[2]][cells.to.keep,]
  gInfoTotals=p1+p2 #phasedCounts[[1]][cells.to.keep,]+phasedCounts[[2]][cells.to.keep,]
  informativeCellsPerTranscript=colSums(gInfoTotals>0)
  genes.to.test=names(informativeCellsPerTranscript)[informativeCellsPerTranscript>informativeCellThresh]
 
  #ef=rnorm(nrow(p2))
  #TO DO: build this (ef) into the ase.Data data structure 
  if(  length(unique(ase.Data$experiment.factor)) > 1 ) { 
      ef=as.factor(ase.Data$experiment.factorNULL[cells.to.keep,])
  } else {
      ef=NULL
  }
  bbin=mclapply(genes.to.test,
    function(g,...){
        ex=cbind(p2[,g],p1[,g])
        print(g)
        if(is.null(ef)){
            m=(glmmTMB(ex~1,family=glmmTMB::betabinomial(link='logit'))) 
        } else{
            m=(glmmTMB(ex~ef,family=glmmTMB::betabinomial(link='logit'))) 
        }
        result=list()
        result$stats=NULL
        result$coefs=NULL
        result$gene = g
        if(!is.na(logLik(m))) {
            result$stats=glance(m)
            result$coefs= tidy(m, effects='fixed', component="cond", exponentiate=F, conf.int=F)
        }
        return(result)
    },p1,p2,ef,mc.cores=threads)
  names(bbin)= sapply(bbin, function(x) x$gene) #genes.to.test[1:240]

  bbin.model.stats=rbindlist(lapply(bbin, function(x) x$stats),idcol='gene', fill=T)
  bbin.model.coefs=rbindlist(lapply(bbin[bbin.model.stats$gene[!is.na(bbin.model.stats$logLik)]], 
                                    function(x) { x$coefs}),
                                        #tidy(x, effects='fixed', component="cond", exponentiate=T, conf.int=T)),
                             idcol='gene', fill=T)
  bbin.model.results=left_join(bbin.model.coefs, bbin.model.stats)
  return(bbin.model.results)
}


#come on broom.mixed, why no function to extract tidy info for disp part of glm fit
tidyDisp=function(m){
    ctable=(coef(summary(m))$disp)
    colnames(ctable)=c('estimate', 'std.error', 'statistic', 'p.value')
    ctable=tibble(data.frame(effect='fixed', component='disp', term=rownames(ctable), 
                      ctable, conf.low=NA, conf.high=NA, stringsAsFactors=F))
    return(ctable)
}


doNbinTest=function(dip,phasedCounts,dip.Assignments,ase.Data, 
                        reduced.gene.set=NULL,
                        nUMI.thresh=10000,informativeCellThresh=64,threads=48) {

  classified.cells=dip.Assignments$diploid_name==dip
  selected.barcodes=dip.Assignments$barcode[dip.Assignments$diploid_name==dip]
  cells.to.keep=ase.Data$numi[classified.cells]<nUMI.thresh 
  p1=phasedCounts[[1]][cells.to.keep,]
  p2=phasedCounts[[2]][cells.to.keep,]
  gInfoTotals=p1+p2 #phasedCounts[[1]][cells.to.keep,]+phasedCounts[[2]][cells.to.keep,]
  informativeCellsPerTranscript=colSums(gInfoTotals>0)
  genes.to.test=names(informativeCellsPerTranscript)[informativeCellsPerTranscript>informativeCellThresh]

  if(  length(unique(ase.Data$experiment.factor)) > 1 ) { 
      ef=as.factor(ase.Data$experiment.factorNULL[cells.to.keep,])
  } else {
      ef=NULL
  }
  ncells=sum(cells.to.keep)
  cellID=factor(as.character(seq(1:ncells))) #length(r)
  setup.vars=list(
        cellIDf=factor(c(cellID,cellID)),
        geno=factor(c(rep('A',ncells),rep('B',ncells))),
        of2=c(log(ase.Data$numi[classified.cells][cells.to.keep]),
               log(ase.Data$numi[classified.cells][cells.to.keep])),
        ef=ef
        )


  nbin=mclapply(genes.to.test[1:48], 
   function(g, ... ){
        print(g)
        #r=rMat[cells.to.keep,g]
        #a=aMat[cells.to.keep,g]
        #efs=experiment.factor[paroi][numi[paroi]<5000]
        #y=c(r,a)#[numi
        y=c(p2[,g],p1[,g])
        cellIDf=setup.vars$cellIDf
        of2=setup.vars$of2
        geno=setup.vars$geno
        ef=setup.vars$ef
        if(is.null(ef)){
            m=glmmTMB(y~geno+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
        }
        else {
            m=glmmTMB(y~ef+geno+(1|cellIDf)+offset(of2), dispformula=~geno, family=nbinom2(link='log'))
        }
        result=list()
        result$stats=NULL
        result$coefs.cond=NULL
        result$coefs.disp=NULL
        result$gene = g
        if(!is.na(logLik(m))) {
            result$stats=glance(m)
            result$coefs.cond= tidy(m, effects='fixed', component="cond", exponentiate=F, conf.int=T)
            result$coefs.disp= tidyDisp(m)

        }
        return(result)
        #if(!is.na(logLik(nbfit))){     nbin[[g]]=summary(nbfit)  } else {nbin[[g]]=NULL }
   },
  p1,p2,setup.vars, mc.cores=48)

 names(nbin)= sapply(nbin, function(x) x$gene) #genes.to.test[1:240]

 nbin.model.stats=rbindlist(lapply(nbin, function(x) x$stats),idcol='gene', fill=T)
 nbin.model.coefs=rbindlist(lapply(nbin[nbin.model.stats$gene[!is.na(nbin.model.stats$logLik)]], 
                                    function(x) { rbind(x$coefs.cond, x$coefs.disp)}),
                             idcol='gene', fill=T)
 nbin.model.results=left_join(nbin.model.stats, nbin.model.coefs, by='gene')
 return(nbin.model.results)
}
