# Build r/qtl2 cross object
# takes genotype calls and genetic map as input
# for now either expecting f2 like data (multi-generation intercross)
# or haploid data (for now using r/qtl 'riself') but should update to new r/qtl2 class that handles haploids?
buildCrossObject=function(graw, gmap.s, cross.type='f2' ) { 
    #chr
    chr=tstrsplit(rownames(graw), ':')[[1]]
    #cell name
    id=colnames(graw)

    # Build r/qtl (v1) cross object to hack input to HMM ------------------
    # hack a r/qtl cross object without the annoying formatting and i/o 
    # r/qtl can handle partially informative genotypes, so let's start there
    # https://www.biostat.wisc.edu/~kbroman/teaching/misc/SISG/2005/hmm_imp_handout.pdf
    # g.recode is markers by segregants
    # g.recode expects '-', 1 or 2 
    # 'C' = notA (het or B)
    # 'D' = notB (A or het)
    # R/qtl codes C and D as 4 and 5 (ok, whatever)
    if(cross.type=='f2'){
        g.split=lapply(split(data.frame(graw+4), chr)  , data.matrix)
    }
    if(cross.type=='riself'){
        g.split=lapply(split(data.frame(graw+1), chr)  , data.matrix)
    }
    #g.split=lapply(split(data.frame(graw+1), chr)  , data.matrix)
    g.split=lapply(g.split, function(x) list(data=t(x)) )
    # all autosomes for now, X chromosome needs to be handled more carefully
    g.split=lapply(g.split, function(x) { class(x)='A'; return(x); })
    # in theory phenotype data could be added here 
    pmat=data.frame(id=id)
    cross=list(geno=g.split, pheno=pmat)
    class(cross)=c(cross.type, 'cross')

    # add genetic map -------------------------------------------------------
    for(uc in names(cross$geno)) {
        cross$geno[[uc]]$map=gmap.s[[uc]]
    }
    #-----------------------------------------------------------------------
    # now turn it into an r/qtl2 object
    cross=convert2cross2(cross)
    return(cross)
}

# should probably just start with this 
split_gmap=function(graw,gmap) {
    chr=tstrsplit(rownames(graw), ':')[[1]]
    names(gmap)=rownames(graw)
    gmap.s=split(gmap, chr)
    return(gmap.s)
}

# given sparseness of genotyping data, consider downsampling markers (e.g. 1 every 5-10 cm for further speedup)
reduce_markersJB=function(gmap.s, gps, min_distance=1){
 
   red_markers=reduce_markers(gmap.s, min_distance=min_distance)

    # come on, this is annoying
    gps.s=mapply(function(x,y){
       nn=match(names(y), dimnames(x)[[3]])
       nn2=dimnames(x)[[3]][nn]
       z=x[,,nn]
       dimnames(z)[[1]]=dimnames(x)[[1]]
       dimnames(z)[[2]]=dimnames(x)[[2]]
       dimnames(z)[[3]]=nn2
       return(z)
    }, x=gps, y=red_markers, SIMPLIFY=F, USE.NAMES=T)
    attr(gps.s, 'crosstype')='f2'
    attr(gps.s,'is_x_chr')=rep(F,6)
    names(attr(gps.s,'is_x_chr'))=as.character(as.roman(1:6))
    attr(gps.s,'alleleprobs')=F
    attr(gps.s,'class')=c('calc_genoprob', 'list')
    return(gps.s)
}

# calculate a FDR threshold as in Smith et al. 2008
calcFDRthresh=function(gscan, gscanPerm,FDR.thresh=c(.20,.10,.05)) {
    max.obsLOD=na.omit(colMaxs(gscan, value=T))
    maxPerms=sapply(gscanPerm, function(x) colMaxs(x, value=T))
    obsPcnt = sapply(seq(1.5, 9, .05), function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = seq(1.5, 9, .05)
    if(sum(obsPcnt)<3) {break}
    # expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(seq(1.5, 9, .05),  
                                 function(thresh) { 
                                        #print(thresh); 
                                        mean(apply(maxPerms, 2, function(ll) {sum(ll>thresh) }) )
                                    } )
    names(expPcnt) = seq(1.5, 9, .05)
    pFDR = expPcnt/obsPcnt
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    pFDR = rev(cummax(rev(pFDR)))
    fdrFX=approxfun(pFDR, seq(1.5,9,.05))
    thresh=fdrFX(FDR.thresh)
    names(thresh)=FDR.thresh
    return(thresh)
}
# report on QTL peaks using r/qtl functions
getPeaks=function(gscan, gmap.s, geneAnnot,threshold){
    
    contig.lengths=sapply(split(geneAnnot$start_position,geneAnnot$chr), max)[-5]
    gcoord.key=cumsum(c(0, contig.lengths))[-7]
    names(gcoord.key)=c(as.character(as.roman(1:5)), 'X')

    peaks=find_peaks(gscan, gmap.s, threshold, drop=1.5, sort_by='lod')
    peaks$peak_pmap=find_marker(gmap.s, peaks$chr, peaks$pos)
    peaks$peak_pmap_L=find_marker(gmap.s, peaks$chr, peaks$ci_lo)
    peaks$peak_pmap_R=find_marker(gmap.s, peaks$chr, peaks$ci_hi)
    peaks=merge(peaks, geneAnnot, by.x='lodcolumn', by.y='external_gene_name', all.x=T, sort=F)
    peaks=peaks[,c(1,11:16,3:10)]
    peaks$peak_pmap_pos=tstrsplit(peaks$peak_pmap, ':',type.convert=T)[[2]]
    peaks$is.local=F
    peaks$is.local[which(peaks$chromosome_name==peaks$chr)]=T
    peaks=peaks[peaks$chromosome_name %in% c(as.character(as.roman(1:5)), 'X'),]
    peaks$chromosome_name=as.factor(peaks$chromosome_name)
    peaks$peak.gcoord=gcoord.key[peaks$chr]+peaks$peak_pmap_pos   
    peaks$t.gcoord=gcoord.key[peaks$chromosome_name]+peaks$start_position   
    return(peaks)
}


