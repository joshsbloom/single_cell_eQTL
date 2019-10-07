getGenoInformativeCounts=function(data.dir, cranger.dir,gmap,parents,log2diff=3){
    GT.cols=grep('GT', names(parents))
    parents$m2=paste0(parents$chr, '_',parents$pos)


    aC=readMM(paste0(data.dir, 'alt_counts.mtx'))
    rC=readMM(paste0(data.dir, 'ref_counts.mtx'))

    # 0 index from vartrix
    variants=read.csv(paste0(data.dir,'out_var.txt'),sep='\t',header=F,stringsAsFactors=F)[,1]
    variants= tstrsplit(variants,'_', type.convert=T)
    variants[[2]]=variants[[2]]+1

    # marker for which we have no corresponding information in the genetic map
    bad.markers=which(!paste0(variants[[1]], '_', variants[[2]]) %in% paste0(gmap$chrom, '_', gmap$ppos))
    variants[[1]]=variants[[1]][-bad.markers]
    variants[[2]]=variants[[2]][-bad.markers]

    ref.counts=rC[-bad.markers,]
    alt.counts=aC[-bad.markers,]
   
    rownames(ref.counts)=paste0(variants[[1]],'_', variants[[2]])
    colnames(ref.counts)=barcodes

    rownames(alt.counts)=rownames(ref.counts)
    colnames(alt.counts)=colnames(ref.counts)

    
    #for haploids for now code crosstype="riself", genotypes as 1 or 2 (missing is 0)
    # setup r/qtl object

    #ugly fix this    
    more.bad.markers=match(rownames(ref.counts) ,parents$m2)
    if(sum(is.na(more.bad.markers))>0){
        ref.counts=ref.counts[-which(is.na(more.bad.markers)),]
        alt.counts=alt.counts[-which(is.na(more.bad.markers)),]
    }

    i1=which(ref.counts>0, arr.ind=T)
    i2=which(alt.counts>0, arr.ind=T)
    i3=rbind(i1,i2)
    informative.variants=unique(sort(i3[,1]))
    informative.cells=unique(sort(i3[,2]))

    #most robust option is to flip ref/alt here!
    # at this step, rQTL.coded needs to be flipped 
    #seg.mat=rQTL.coded-1
   

    #sanity check ref/alt variants here 
    #grep(',', parents$alt)
    
    # the match function here could be problematic if vcf contains multiple entries at a given position 
    # double check this step
    parents.s=parents[match(rownames(ref.counts), parents$m2),] #parent
    multivar=(grepl(',', parents.s$alt) | grepl(',', parents.s$ref))
    
    parent.phasing=parents.s[,GT.cols]
    p1.ref=parent.phasing[,1]<parent.phasing[,2]
  
    #double check this 
    p1=rbind(ref.counts[p1.ref,],alt.counts[!p1.ref,])
    p1=p1[rownames(ref.counts),]
   
    p2=rbind(ref.counts[!p1.ref,],alt.counts[p1.ref,])
    p2=p2[rownames(ref.counts),]

    # double check this 
    p1[which(multivar),]=0
    p2[which(multivar),]=0
    #ref.counts
    #nn=2761 #2993 #4173
    #plot(p1[,nn], ylim=c(-10,10), type='h')
    #points(-p2[,nn], ylim=c(-10,10), type='h')
    
    ref.counts=p1
    alt.counts=p2
    
    #distortions could be AF/ASE/false positive variants/
    distorted=cbind(log2(rowSums(ref.counts)+1), log2(rowSums(alt.counts)+1))
    filter.distorted=(abs(distorted[,1]-distorted[,2])>log2diff)
    plot(distorted[,1], distorted[,2], col=1+filter.distorted, sub=sum(filter.distorted)/nrow(distorted))
    # I think we can afford to be aggressive with filtering here 
    
    ref.counts[which(filter.distorted),]=0
    alt.counts[which(filter.distorted),]=0

    return(list(ref.counts=ref.counts,alt.counts=alt.counts))
}
