getGenoInformativeCounts=function(vartix.alt, vartrix.ref, gmap,
                                      vartrix.parental.ref.file,
                                      vartrix.parental.alt.file,
                                      vartrix.parental.monocle.file    ){
    # counts for alternative allele
    aC=readRDS(vartrix.alt)
    # counts for reference allele
    rC=readRDS(vartrix.ref)

    # read in vartrix output variants 
    variants=rownames(aC) 

    #always some bullshit, Eyal, you are missing 2 markers from your imputed map?
    bad.markers=which(!variants %in% gmap$marker)

    variants=variants[-bad.markers]
    N2.counts=rC[-bad.markers,]
    CB.counts=aC[-bad.markers,]
    
    i1=which(N2.counts>0, arr.ind=T)
    i2=which(CB.counts>0, arr.ind=T)
    i3=rbind(i1,i2)
    informative.variants=unique(sort(i3[,1]))

    #logic here is that ~1/3 of sites have some information, let's just run HMM etc on those sites 
    #with sparse matrices, indexing doesn't change ordering, be careful though
    N2.counts=N2.counts[informative.variants,]
    CB.counts=CB.counts[informative.variants,]
    #originally stopped here 
    ##########################################################


    #####additional fitler on parental vartrix data 
    aC=readRDS(vartrix.alt)
    refC=readRDS(vartrix.parental.ref.file)
    rownames(refC)=rownames(aC)
    altC=readRDS(vartrix.parental.alt.file) 
    rownames(altC)=rownames(aC)
    refC=refC[rownames(N2.counts),]
    altC=altC[rownames(N2.counts),]

    par.monocle=readRDS(vartrix.parental.monocle.file)
    par.monocle.df=colData(par.monocle)@listData

    refC=refC[,!is.na(par.monocle.df$Strain)]
    altC=altC[,!is.na(par.monocle.df$Strain)]
    par.monocle.df=data.frame(par.monocle.df, stringsAsFactors=F)[!is.na(par.monocle.df$Strain),]
    N2.ref=refC[,par.monocle.df$Strain=='N2']
    N2.alt=altC[,par.monocle.df$Strain=='N2']

    N2.rc=rowSums(N2.ref)
    N2.ac=rowSums(N2.alt)
    bad.N2=(log2(N2.ac+1)>log2(N2.rc+1)+1.5)
    par(mfrow=c(1,3))
    plot(log2(N2.rc), log2(N2.ac), col=bad.N2+1, main='cells classified as N2', xlab='N2 counts', ylab='CB counts')

    CB.ref=refC[,par.monocle.df$Strain=='CB4856']
    CB.alt=altC[,par.monocle.df$Strain=='CB4856']
    CB.rc=rowSums(CB.ref)
    CB.ac=rowSums(CB.alt)
    bad.CB=(log2(CB.rc+1)>log2(CB.ac+1)+1.5)
    plot(log2(CB.ac), log2(CB.rc), col=bad.CB+1, main='cells classified as CB', xlab='CB counts', ylab='N2 counts')

    N2.rc.n=rowSums(N2.ref/(colSums(N2.ref)+colSums(N2.alt)))
    N2.ac.n=rowSums(N2.alt/(colSums(N2.ref)+colSums(N2.alt)))

    CB.rc.n=rowSums(CB.ref/(colSums(CB.ref)+colSums(CB.alt)))
    CB.ac.n=rowSums(CB.alt/(colSums(CB.ref)+colSums(CB.alt)))

    offn=.5
    diff.exp.par=abs(log2(N2.rc.n+offn)- log2(CB.ac.n+offn))>3
    more.bad.variants=diff.exp.par | bad.N2 | bad.CB

    plot(log2(N2.rc.n+offn), log2(CB.ac.n+offn), xlab='N2', ylab='CB',col=more.bad.variants+1, main='comparing parental counts')
    abline(0,1)
    abline(lm(log2(CB.ac.n+offn)~log2(N2.rc.n+offn)), col='red')
    lm(log2(CB.ac.n+offn)~log2(N2.rc.n+offn))

    N2.counts=N2.counts[-which(more.bad.variants),]
    CB.counts=CB.counts[-which(more.bad.variants),]
    return(list(N2.counts=N2.counts, CB.counts=CB.counts))
} 
