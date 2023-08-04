 build.gcoord.key =function(filein) {
    ref.seq=read.fasta(filein)
    contig.lengths=sapply(ref.seq, length)
    names(contig.lengths)=paste0('chr', c(as.character(as.roman(1:16)), 'Mito'))
    gcoord.key=cumsum(c(0, as.vector(contig.lengths)))[-18]
    names(gcoord.key)=names(contig.lengths)
    return(gcoord.key)
}

extractGeneticMapDf=function(cross){
    preliminary.map=sapply(cross$geno, function(x) x$map)
    # we need to add some jitter here 
    gmap=data.frame(marker=as.vector(unlist(sapply(preliminary.map, names))),
                    gpos=as.vector(unlist(sapply(preliminary.map, function(x)x))),stringsAsFactors=F)
    gmap2=tstrsplit(gmap$marker,'_')
    gmap$chrom=gmap2[[1]]
    gmap$ppos=as.numeric(gmap2[[2]])
    return(gmap) 
} 


#defunct 
#additional helper function
recode.as.parental=function(seg.mat, parents) { 
    GT.cols=grep('GT', names(parents))
    z=rep(NA,nrow(seg.mat))
    names(z)=rownames(seg.mat)
    g.recode=apply(seg.mat, 2, function(x) {
        z[which(x==parents[,GT.cols[1]])]=1
        z[which(x==parents[,GT.cols[2]])]=2
        return(z)
    })
    return(g.recode)
}

#convert marker positions into granges object
getMarkerGRanges=function(g.counts) {

    g2m=tstrsplit(rownames(g.counts[[1]]),'_')
    markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*"))
    markerGR$name=rownames(g.counts[[1]])
    return(markerGR)
}

getCisMarker=function(sgd.genes, m.granges, counts) {
    Yre=t(counts)
    markerGR=m.granges
    sgd.genes=sgd.genes[match(colnames(Yre), sgd.genes$Name),]
    tag=sgd.genes

    cisMarkers=data.frame(transcript=sgd.genes$Name, marker=m.granges$name[nearest(tag,markerGR)],stringsAsFactors=F)
    cisMarkers=na.omit(cisMarkers)
    return(cisMarkers)
    #indices of cis eQTL in 
    #rows are transcripts, columns are markers
    #cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
}



getASEcounts=function(g.counts, sgd.genes, m.granges, threads=16) {
    # can split this back out if we want to keep track of counts from individual variants     
    fO=findOverlaps(sgd.genes, m.granges)
    mgsplit=split(m.granges$name[subjectHits(fO)], sgd.genes$Name[queryHits(fO)])
    print('calc counts per gene, allele 1')
    # could restructure these as sparse matrices if necessary
    ref.ASEcounts=mcmapply( function(x) {
            if(length(x)>1) {
                colSums(g.counts[[1]][x,])
            } else {
                g.counts[[1]][x,]
            }
            #colSums(g.counts[[2]][x,]) 
       },mgsplit,mc.cores=threads)

    print('calc counts per gene, allele 2')
    alt.ASEcounts=mcmapply( function(x) {
            if(length(x)>1) {
                colSums(g.counts[[2]][x,])
            } else {
                g.counts[[2]][x,]
            }
            #colSums(g.counts[[2]][x,]) 
       }, mgsplit,mc.cores=threads)
        return(list(ref=ref.ASEcounts, alt=alt.ASEcounts))
}

getCCinfo=function(cell.cycle.classification.file,cell.cycle.annot.file) {
    cc.cl=read.delim(cell.cycle.classification.file, header=F, sep=' ', stringsAsFactors=F)
    names(cc.cl)=c('barcode', 'nid')
    cc.cl$cid=NA
    cc.annot=read.delim(cell.cycle.annot.file, header=F, sep=':', stringsAsFactors=F)
    cc.annot.list=strsplit(cc.annot[,2], ',')
    names(cc.annot.list)=cc.annot[,1]
    cc.annot1=cbind(rep(names(cc.annot.list), sapply(cc.annot.list, length)), as.vector(unlist(cc.annot.list)))

    cc.cl$cid=cc.annot1[,1][match(as.character(cc.cl$nid), cc.annot1[,2])]
    #and order the factor 
    cc.cl$cid=factor(cc.cl$cid, levels=cc.annot[,1])
    return(cc.cl)
}

# type can be viterbi or genoprobs
getSavedGenos=function(chroms, results.dir, type='viterbi'){
         if(type=='viterbi') {
             fname='_viterbi.RDS' 
             vg=do.call('cbind', lapply(chroms, function(coii) {
             v.file=paste0(results.dir, coii, fname)
             vit= readRDS(v.file) } ))
             return(vg) 
         }

         if(type=='genoprobs') {
             fname='.RDS'
             vg=do.call('cbind', lapply(chroms, function(coii) {
             v.file=paste0(results.dir, coii, fname)
             vit= readRDS(v.file)[,2,]
             return(vit)
             })) 
             return(vg)
         }
}


getFDRfx=function(rff, Yfe,Gf,nperm=5){
     
    max.obsLOD=apply(abs(rff),1,max)
    vint=seq(0.001,max(max.obsLOD)+.001, .001)
    # no point keeping all this, just retain max stats 
    ll=replicate(nperm, {
                    nind=sample(1:nrow(Yfe))
                    rowMaxs(abs(crossprod(Yfe[nind,],Gf)/(nrow(Yfe)-1)),value=T)
    })

    #ll=apply(abs(permR),3, function(x) apply(x,1,max))
    #rm(permR)

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
    pFDR[1]=pFDR[1]-(1e-3)
    fdrFX=approxfun(pFDR, vint, ties=mean)
    rtoFDRfx=approxfun(vint,pFDR, ties=mean)

    return(list(fdrFX=fdrFX,rtoFDRfx=rtoFDRfx))
}

get1Dpeaks=function(chroms, r, LOD, wgFDR, m.granges, sgd.genes) {
    cPeaks=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(r))
        rS=r[,moi]
        #mstat=apply(rS,1,max)
        mstat.pos=apply(abs(rS),1,which.max)
        lookup=cbind(1:length(mstat.pos), mstat.pos)
        mstat=rS[lookup]
        LODstat=LOD[,moi][lookup]
        CIs=matrix(NA,length(mstat),2)
        cmoi=colnames(LOD[,moi])
        for(peak in 1:length(mstat)){
            CIs[peak,]=cmoi[range(which(LOD[peak,moi]>LODstat[peak]-1.5))]
        }

        cPeaks[[cc]]=data.frame( transcript=rownames(rS),
                                 peak.marker=colnames(rS)[mstat.pos],
                                 CI.l=CIs[,1],
                                 CI.r=CIs[,2],
                                 LOD=LODstat, 
                                 sbeta=mstat,
                                 FDR=wgFDR[[2]](abs(mstat)) )
        cPeaks[[cc]]$FDR[is.na(cPeaks[[cc]]$FDR)]=1
    }
    cPeaks=rbindlist(cPeaks, idcol='chr')
    cPeaks$mchr=as.character(seqnames(m.granges))[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mpos=(start(m.granges))[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mgcoord =m.granges$gcoord[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mgcoordL=m.granges$gcoord[match(cPeaks$CI.l, m.granges$sname)]
    cPeaks$mgcoordR=m.granges$gcoord[match(cPeaks$CI.r, m.granges$sname)]
    cPeaks$mposL   =(start(m.granges))[match(cPeaks$CI.l, m.granges$sname)]
    cPeaks$mposR   =(start(m.granges))[match(cPeaks$CI.r, m.granges$sname)]

    cPeaks$tchr=as.character(seqnames(sgd.genes))[match(cPeaks$transcript, sgd.genes$Name)]
    cPeaks$tpos=start(sgd.genes)[match(cPeaks$transcript, sgd.genes$Name)]
    cPeaks$tgcoord=sgd.genes$gcoord[match(cPeaks$transcript, sgd.genes$Name)]
    return(cPeaks)
}


mvn.scanone=function(g.s, Y,add.cov=NULL, roi1=NULL) {
           ldetRSSf.1=rep(NA, ncol(g.s))
           if(!is.null(roi1)) { ss=roi1 } else { ss=1:ncol(g.s) } 
           for(i in min(ss):max(ss)) {
               #ldetRSSf.1[i]=det_AtA(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               if(is.null(add.cov) ) {
                RSSf.1=crossprod(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               } else {
                RSSf.1=crossprod(residuals(.lm.fit(cbind(add.cov, g.s[,i]), Y)))
               }
               ldetRSSf.1[i]=determinant(RSSf.1, logarithm=T)$modulus
           }
           return(ldetRSSf.1)
}




getMVpeaks=function(chroms, r, LOD, mLOD,wgFDR,m.granges,sgd.genes) {
    cPeaks=list()
    mPeaks=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(r))
        rS=r[,moi]
        rMOI=mLOD[moi]
        #mstat=apply(rS,1,max)
        mstat.pos=which.max(rMOI) #),1,which.max)
        lookup=cbind(1:nrow(rS), mstat.pos)
        mstat=rS[lookup]
        LODstat=LOD[,moi][lookup]
        mLODstat=mLOD[moi][mstat.pos]

        CIs=matrix(NA,length(mstat),2)
        cmoi=colnames(LOD[,moi])
        for(peak in 1:length(mstat)){
            CIs[peak,]=cmoi[range(which(LOD[peak,moi]>LODstat[peak]-1.5))]
        }
        mci=cmoi[range(which(mLOD[moi]>mLODstat-1.5))]
        mPeaks[[cc]]=data.frame(mLOD=mLODstat,
                                peak.marker=colnames(rS)[mstat.pos],
                                 CI.l=mci[1],
                                 CI.r=mci[2])
        cPeaks[[cc]]=data.frame( transcript=rownames(rS),
                                 peak.marker=colnames(rS)[mstat.pos],
                                 CI.l=CIs[,1],
                                 CI.r=CIs[,2],
                                 LOD=LODstat, 
                                 sbeta=mstat,
                                 FDR=wgFDR[[2]](abs(mstat)) )
        cPeaks[[cc]]$FDR[is.na(cPeaks[[cc]]$FDR)]=1
    }
    cPeaks=rbindlist(cPeaks, idcol='chr')

    cPeaks$mchr=as.character(seqnames(m.granges))[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mpos=(start(m.granges))[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mgcoord=m.granges$gcoord[match(cPeaks$peak.marker, m.granges$sname)]
    cPeaks$mgcoordL=m.granges$gcoord[match(cPeaks$CI.l, m.granges$sname)]
    cPeaks$mgcoordR=m.granges$gcoord[match(cPeaks$CI.r, m.granges$sname)]
    cPeaks$mposL   =(start(m.granges))[match(cPeaks$CI.l, m.granges$sname)]
    cPeaks$mposR   =(start(m.granges))[match(cPeaks$CI.r, m.granges$sname)]


    cPeaks$tchr=as.character(seqnames(sgd.genes))[match(cPeaks$transcript, sgd.genes$Name)]
    cPeaks$tpos=start(sgd.genes)[match(cPeaks$transcript, sgd.genes$Name)]
    cPeaks$tgcoord=sgd.genes$gcoord[match(cPeaks$transcript, sgd.genes$Name)]

    mPeaks=rbindlist(mPeaks, idcol='chr')
    mPeaks$mpos=(start(m.granges))[match(mPeaks$peak.marker, m.granges$sname)]
    mPeaks$mgcoord=m.granges$gcoord[match(mPeaks$peak.marker, m.granges$sname)]
    mPeaks$mgcoordL=m.granges$gcoord[match(mPeaks$CI.l, m.granges$sname)]
    mPeaks$mgcoordR=m.granges$gcoord[match(mPeaks$CI.r, m.granges$sname)]
    mPeaks$mposL   =(start(m.granges))[match(mPeaks$CI.l, m.granges$sname)]
    mPeaks$mposR   =(start(m.granges))[match(mPeaks$CI.r, m.granges$sname)]

    return(list(cPeaks=cPeaks, mPeaks=mPeaks))
}

plot2D=function(peaks1D, FDR.thresh=.05, gcoord.key=gcoord.key, experiment=experiment) {
    plot.thresh=FDR.thresh
    par(xaxs='i', yaxs='i')
    plot(peaks1D$mgcoord[peaks1D$FDR<plot.thresh], peaks1D$tgcoord[peaks1D$FDR<plot.thresh],
         ylab='transcript position',
         xlab='QTL position', main=experiment, xaxt='n', yaxt='n', pch=20, 
        sub=paste('FDR<', FDR.thresh*100, '%'))
    axis(1, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    axis(2, at=rowMeans(cbind(gcoord.key[-1], gcoord.key[-17])), labels=names(gcoord.key)[-17])
    abline(v=gcoord.key, lty=2, col='grey')
    abline(h=gcoord.key, lty=2, col='grey')
}




# parallelized standardise function
standardise2 <- function (x, center = TRUE, scale = TRUE) {
  if ( center & scale ) {
    y <- t(x) - Rfast::colmeans(x,parallel=T)
    y <- y / sqrt( Rfast::rowsums(y^2) ) * sqrt( (dim(x)[1] - 1) )
	y <- t(y)
  } else if ( center & !scale ) {
    m <- Rfast::colmeans(x,parallel=T)
    y <- eachrow(x, m, oper ="-" )
  } else if ( !center & scale ) {
    s <- Rfast::colVars(x, std = TRUE,parallel=T)
    y <- eachrow(x, s, oper = "/")
  }
  y
} 

