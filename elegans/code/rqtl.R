#adopted from r/qtl
jitterGmap=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- nrow(themap[[i]])
         themap[[i]]$map <- themap[[i]]$map + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}


# hard call genotypes and encode for r/qtl hmm
encodeForRQTL_F2=function(N2.counts, CB.counts) { 
    # Code input to r/qtl
    # see hmm_f2.c from r/qtl
    #  be careful here 
    #               1    2    3    4    5
    #genotypes = c("A", "H", "B", "C", "D")
    # see papers from Lynch and Maruki if we want to extract a bit more infomation here
    # N2 as AA
    # CB as BB

    rQTL.coded=matrix(NA, nrow(N2.counts), ncol(N2.counts))
    rownames(rQTL.coded)=rownames(N2.counts)
    colnames(rQTL.coded)=colnames(N2.counts)
    mode(rQTL.coded)="integer"
    # a bit of memory reduction by changing mode

    # having some issues doing logic operations on R sparse matrices
    # here is the workaround 
    N2=matrix(0, nrow(N2.counts), ncol(N2.counts))
    rownames(N2)=rownames(N2.counts)
    colnames(N2)=colnames(N2.counts)
    mode(N2)='integer'
    CB=N2


    N2gt0=which(N2.counts>0)
    CBgt0=which(CB.counts>0)

    N2[N2gt0]=as.vector(N2.counts[N2gt0])
    CB[CBgt0]=as.vector(CB.counts[CBgt0])

    hets= N2gt0[N2gt0 %in% CBgt0]
    rQTL.coded[hets]=2

    # C is not AA (B or het)
    CC=CBgt0[!(CBgt0 %in% N2gt0)]
    rQTL.coded[CC]=5

    # D is not BB (A or het)
    DD=N2gt0[!(N2gt0 %in% CBgt0)] 
    rQTL.coded[DD]=4

    n2.10=which(N2.counts>10)
    cb.10=which(CB.counts>10)

    A=n2.10[!(n2.10 %in% CBgt0)]
    B=cb.10[!(cb.10 %in% N2gt0)]
    rQTL.coded[A]=1
    rQTL.coded[B]=3

    #re-evaluate het calls
    N2h=N2[hets]
    CBh=CB[hets]

    # n = total counts
    # k = ref counts
    # r = prior that call is het
    # e = genotyping error rate
    #pHet choose(n,k)*(1/2^n)*r

    n=N2h+CBh
    k=N2h
    r=.5
    ee=.0001
    #n=2
    #k=0

    #for example
    pHet = choose(n,k)*(1/(2^n))*   r
    pAA  = dbinom(k,n,ee)       *   (1-r)/2
    pRR  = dbinom(n-k,n,ee)     *   (1-r)/2

    reev_het=apply(cbind(pHet,pAA,pRR),1,which.max)
    reev_het[is.na(pHet)]=1
    rQTL.coded[hets[which(reev_het==2)]]=3
    rQTL.coded[hets[which(reev_het==3)]]=1

   # tcount=matrix(0, nrow(rQTL.coded), 5)
   # for(i in 1:nrow(rQTL.coded)){
   #     if(i%%1000==0) {print(i)}
   #     gt=(table(rQTL.coded[i,]))
   #     tcount[i,c(as.numeric(names(gt)))]=as.vector(gt)
   # }

   #3 CC    5    NC/CC 
   #1 NN    4 NN/NC
   #2 NC
    return(rQTL.coded)
}

#turn into r/QTL cross object  and then rQTL2 cross2 object
buildCross2=function(rQTL.coded, gmap.s, counts) {
    #get this straightened out!
    # r/qtl hack formatting bullshit slow as fuck
    g.split=lapply(split(data.frame(rQTL.coded),gmap.subset$chrom), data.matrix)
    g.split=lapply(g.split, function(x) list(data=t(x)) )
    g.split=lapply(g.split, function(x) { class(x)='A'; return(x); })
    gc()
    pmat=data.frame(id=colnames(counts))
    cross=list(geno=g.split, pheno=pmat)
    class(cross)=c('f2', 'cross')

    # add genetic map -------------------------------------------------------
    for(uc in names(cross$geno)) {
            cross$geno[[uc]]$map=gmap.s[[uc]]$map
            names(cross$geno[[uc]]$map)=paste0(gmap.s[[uc]]$chrom,'_', gmap.s[[uc]]$pos)
    }
    #-----------------------------------------------------------------------
    rm(rQTL.coded)
    rm(g.split)
    #cross=jittermap(cross)
    cross2=convert2cross2(cross)
    return(cross2)
    #saveRDS(cross2, file='/data/single_cell_eQTL/elegans/XQTL_F4_2/cross2_v2.RData')
    #cross2=readRDS('/data/single_cell_eQTL/elegans/XQTL_F4_2/cross2_v2.RData')
}


##Run r/qtl HMM
rQTLHMM=function(cross2.file,hmm.out.dir) {
    #cross2.file= '/data/single_cell_eQTL/elegans/XQTL_F4_2/cross2_v2.RData'                
    #hmm.out.dir='/data/single_cell_eQTL/elegans/XQTL_F4_2/hmm_v2/'
    cross2=readRDS(cross2.file)

    nind=nrow(cross2$cross_info)
    cc=cut(1:nind, 100)

    dir.create(hmm.out.dir)
    #mgeno=list()
    iter=1
    for(ccc in unique(cc)){
        kind=cc==ccc
        print(max(which(kind)))
        c2=cross2[which(kind),]
        gps=calc_genoprob(c2, error_prob=.005,lowmem=F, quiet=F, cores=6)
        saveRDS(gps, file=paste0(hmm.out.dir, iter))
        iter=iter+1
    }
}

# HMM diagnostic plots
    


#    G=readRDS('/data/single_cell_eQTL/elegans/XQTL_F4_2/additive.geno_v2.RDS')
#    
#    hardcall=cross2$geno[[coii]][indiv,]
#    par(oma=c(2,2,2,2))
#    plot(gdist,-N2.counts[mind,indiv], type='h', ylim=c(-5,5),col=hardcall, main =rownames(G)[indiv], ylab='-N2, + CB')
#    abline(h=0)
#    points(gdist,CB.counts[mind,indiv], type='h', ylim=c(-5,5),col=hardcall)
#    legend('topright', c('CC', 'NC', 'NN', 'NN or NC', 'NC or CC'), fill=c('green', 'red', 'black', 'blue', 'cyan'))
#
#    par(new = T)
#    plot(gdist,G[indiv,mind], col='red', type='l', axes=F, ylim=c(0,1), lwd=2, ylab='')
#    points(gdist,(posteriorProb[indiv,3,]+.5*posteriorProb[indiv,2,]), col='purple', type='l', lwd=2)
#    axis(side = 4)
#    mtext(side = 4, line = 3, 'Posterior Prob(CB)')
#    readline()

#black is NN
#red is NC
#green is CC
#blue is NN or NC
#cyan is NC or CC
#------------------------------------------------------------------------------------------------------------------------------





