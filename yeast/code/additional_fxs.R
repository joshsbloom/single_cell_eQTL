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

