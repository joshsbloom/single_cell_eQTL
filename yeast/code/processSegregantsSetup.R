library(monocle)
library(Matrix)
library(rtracklayer)
library(data.table)
library(qtl)
library(qtl2)
library(parallel)
library(Rfast)
library(seqinr)
library(GenomicRanges)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(foreach)
library(doMC)
library(glmmTMB)
library(ggplot2)
library(mgcv)
library(Rfast2)
library(fastglm)
library(MASS)
library(irlba)
library(nebula)
library(dplyr)

# for parallelizing the HMM
ncores=16
registerDoMC(cores=ncores)

code.dir='/data/single_cell_eQTL/yeast/code/'
source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))

source(paste0(code.dir, 'ASE_fxs.R'))

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
  colnames(y)=colnames(x)
  rownames(y)=rownames(x)
  y
} 


fitNBCisModel=function(cMsubset, mmp1, Yk, Gsub, resids=NULL) {
    mmp0=mmp1[,-1]

    #nbR=Yk
    #nbR[is.numeric(nbR)]=NA
    ##nbP=nbR


    cisNB=list() 
    for(r in 1:nrow(cMsubset)){
        print(r)
        gn=as.character(cMsubset$transcript[r])
        p=as.character(cMsubset$marker[r])
        cmodel=cbind(mmp0, Gsub[,p])
        nbr=negbin.reg(Yk[,gn], cmodel,maxiters=500)
        if(!is.na(nbr$info[4]) & (nbr$info[4]<100 & !is.na(nbr$info[3])) ){
             theta.est=nbr$info[4]
             coef.est=nbr$be[,1]
        } else {
            nbr=mgcv::gam(Yk[,gn]~cmodel, family=mgcv::nb , method="ML")
            theta.est=nbr$family$getTheta(T)
            coef.est=as.vector(coef(nbr))
        }
        XX=cbind(mmp1, Gsub[,p])
        msize=ncol(XX)
        fnbrF=fastglm(XX, Yk[,gn], start=coef.est, family=negative.binomial(theta=theta.est,link='log'), maxit=500)
        ##pseudo r^2 statistics 
        ##library(performance)
        ##class(fnbrF)='glm'
        ##r2(fnbrF)
        fnbrN=fastglm(mmp1,Yk[,gn], start=coef.est[-msize],  family=negative.binomial(theta=theta.est,link='log'), maxit=500)
        #if(resids=='deviance'){       nbR[,gn]=residuals(fnbrF,'deviance') }
        #if(resids=='pearson') {       nbR[,gn]=residuals(fnbrF,'pearson') }
        cisNB[[gn]]$transcript=gn
        cisNB[[gn]]$lmarker=p
        cisNB[[gn]]$theta=theta.est
        cisNB[[gn]]$negbin.beta=as.numeric(coef(fnbrF)[msize])
        cisNB[[gn]]$negbin.se=as.numeric(fnbrF$se[msize])
        cisNB[[gn]]$LLik=as.numeric(logLik(fnbrF))
        cisNB[[gn]]$negbin.LRS=as.numeric(-2*(logLik(fnbrN)-logLik(fnbrF)))
        cisNB[[gn]]$negbin.p=pchisq( cisNB[[gn]]$negbin.LRS,1, lower.tail=F) 
         #-2*(nmodel-plr)
        cisNB[[gn]]$fmodelBs=coef(fnbrF)
        
        if(resids=='pearson') {       cisNB[[gn]]$pearsonResiduals=residuals(fnbrF,'pearson') }
        
        #print(cisNB[[gn]])
    }

   # nbR.var=colVars(nbR)
   # nbR=nbR[,-which(is.na(nbR.var))]
    return(cisNB) #, nbR=nbR))
}

# formula for calculating the negative binomial glm log likelihood 
nbLL=function (y, mu, theta)  {
    return( -sum ( (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta + y) ))
}
# fast calculation of negative binomial LL across the genome if theta is known 
domap=function(gn, ss=F,...) { 
       print(gn)
       theta.est=thetas[gn]
       YY=Y[,gn]
       nbn=negative.binomial(theta=theta.est,link='log')
      
       fnbrN=fastglmPure(DM, YY, family=nbn)
       nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)

       LRS=rep(NA,ncol(Gsub))
       names(LRS)=colnames(Gsub)

       if(ss){
            Betas=LRS
            SEs=LRS
       }
       #initialize the last column of DM to be 1s
       XXn=cbind(DM,1)
       idx=1:ncol(Gsub)
       gcovn=(ncol(XXn))
       for(gx in idx) { 
           #colnames(Gsub)){
           XXn[,gcovn]=Gsub[,gx]
           fnbrF=fastglmPure(XXn, YY,  family=nbn)
           LRS[gx]=-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))
           if(ss){
            Betas[gx]=as.numeric(fnbrF$coefficients[gcovn])
            SEs[gx]=fnbrF$se[gcovn]
           }
        }
       if(ss) {
        return(list(LRS=LRS,Betas=Betas,SEs=SEs))
       } else{

       return(LRS)
       }
}

#logistic regression mapping (for discrete cycle classification assignments ) 
domap_logistic=function(gn, ss=F, ...) {
    YY=cc.incidence[,gn]
    DM=mmp1
    fnbrN=fastglmPure(DM, YY, family=binomial())
    #nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)

    LRS=rep(NA,ncol(Gsub))
    names(LRS)=colnames(Gsub)
 
    if(ss){
            Betas=LRS
            SEs=LRS
       }

    #initialize the last column of DM to be 1s
    XXn=cbind(DM,1)
    idx=1:ncol(Gsub)
    gcovn=(ncol(XXn))
    for(gx in idx) { 
       #colnames(Gsub)){
        XXn[,gcovn]=Gsub[,gx]
        fnbrF=fastglmPure(XXn, YY,  family=binomial())
        #can just use deviances here 
        LRS[gx]=fnbrN$deviance-fnbrF$deviance
        if(ss){
            Betas[gx]=as.numeric(fnbrF$coefficients[gcovn])
            SEs[gx]=fnbrF$se[gcovn]
           }
        }
       if(ss) {
        return(list(LRS=LRS,Betas=Betas,SEs=SEs))
       } else{

       return(LRS)
       }
}







           #-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))

#cc=read.csv(paste0(reference.dir, 'cell_cycle.csv')) #read.csv(paste0(

LDprune=function(geno, m.granges, cor.cutoff=.999) {
    G=standardise2(vg)
    #prune markers 
    chr.ind=split(1:length(m.granges), as.character(seqnames(m.granges)))
    G.chr=lapply(chr.ind, function(x) G[,x])
    Gcor.chr=lapply(G.chr, function(x) crossprod(x)/(nrow(x)-1))
    fcm.chr=lapply(Gcor.chr, function(x) caret::findCorrelation(x,cutoff=cor.cutoff, verbose=F) )
    fcm=as.vector(unlist(mapply(function(x,y) x[y],  sapply(G.chr, colnames), fcm.chr)))
    m.to.keep=which(!colnames(vg)%in% fcm) #colnames(vg)[fcm])
    Gsub=vg[,m.to.keep]
    rm(Gcor.chr)
    m.to.keep=match(colnames(Gsub), colnames(vg))
    markerGRr=m.granges[m.to.keep]
    return(list(Gsub=Gsub, markerGRr=markerGRr))
}


#overwrite these functions

getCisMarker=function(sgd.genes, m.granges, counts) {
    Yre=t(counts)
    markerGR=m.granges
    sgd.genes=sgd.genes[match(colnames(Yre), sgd.genes$gene_id),]
    tag=sgd.genes

    cisMarkers=data.frame(transcript=sgd.genes$gene_id, 
                          marker=m.granges$sname[nearest(tag,markerGR)],
                          stringsAsFactors=F)
    cisMarkers=na.omit(cisMarkers)
    return(cisMarkers)
    #indices of cis eQTL in 
    #rows are transcripts, columns are markers
    #cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
}

doCisNebula2=function(x,...) { # mmp1,counts,Gsub){
        mmpN=cbind(mmp1, Gsub[,x$marker[1]])
        yind=match(x$transcript, rownames(counts))
        NBcis=nebula(counts[yind,], mmpN[,1], pred=mmpN)
        
        #NBcis$summary$gene=rownames(counts)[yind]
        return(NBcis)
    }

#for the original ByxRM with repeated measures
doCisNebula3=function(x,...) { # mmp1,counts,Gsub){
       
        mmpN=cbind(mmp1, Gsub[,x$marker[1]])
        yind=match(x$transcript, rownames(counts))
       
        df=data.frame(id=factor(best_match_seg))
      
        mmpN=mmpN[order(df$id),]
        countmat=counts[yind,order(df$id)]
        df.id=df[order(df$id),]

        NBcis=nebula(count=countmat, id=df.id, pred=mmpN) 

        #NBcis=nebula(counts[yind,], mmpN[,1], pred=mmpN)
        NBcis$summary$gene=rownames(counts)[yind]
        return(NBcis)
    }


as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

getGlobalPeaks=function( LODr,markerGR,transcript.data, 
                        chroms=paste0('chr', as.roman(1:16)),
                        fdrfx.fc=NULL,Betas=NULL,SEs=NULL) {
    cPeaksT=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(LODr))
        rS=LODr[,moi]
        #mstat=apply(rS,1,max)
        mstat.pos=Rfast::rowMaxs(abs(rS),value=F) #1,which.max)
        lookup=cbind(1:length(mstat.pos), mstat.pos)
        mstat=rS[lookup]
        LODstat=LODr[,moi][lookup]
        CIs=matrix(NA,length(mstat),2)
        cmoi=colnames(LODr[,moi])
        #this could use a speedup

        for(peak in 1:length(mstat)){
            LODrv=LODr[peak,moi]
            CIs[peak,]=cmoi[range(which(LODrv>LODstat[peak]-1.5))]
        }   
        cPeaksT[[cc]]=data.frame( transcript=rownames(rS),
                             peak.marker=colnames(rS)[mstat.pos],
                             CI.l=CIs[,1],
                             CI.r=CIs[,2],
                             LOD=LODstat)
       
        if(!is.null(Betas)){
            if(is.list(Betas)){
                 for(n in names(Betas)){
                    cPeaksT[[cc]][[paste0('Beta_',n)]]= Betas[[n]][,moi][lookup]
                 }
            }  else {
                cPeaksT[[cc]]$Beta=Betas[,moi][lookup]
            }
        }
        if(!is.null(SEs)){
            if(is.list(SEs)){
                for(n in names(SEs)){
                    cPeaksT[[cc]][[paste0('SE_',n)]]= SEs[[n]][,moi][lookup]
                 }

            } else {
            cPeaksT[[cc]]$SE=SEs[,moi][lookup]
            }
        }

        if(!is.null(fdrfx.fc)){
           cPeaksT[[cc]]$FDR=fdrfx.fc[['g.rtoFDRfx']](LODstat)
           cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
           cPeaksT[[cc]]$FDR[(cPeaksT[[cc]]$FDR)>1]=1
        }
    }
    #inverse variance weights
    #        w=1/SE^2
    #        ivw.mean=sum(w*B)/sum(w)
    #        ivw.se=sqrt(1/sum(w))

    cP=rbindlist(cPeaksT, idcol='chrom')
    cP$peak.marker=as.character(cP$peak.marker)
    cP$chr=tstrsplit(cP$peak.marker, '_')[[1]]
    cP$pos=as.numeric(tstrsplit(cP$peak.marker, '_')[[2]])
    cP$CI.l=as.character(cP$CI.l)
    cP$CI.r=as.character(cP$CI.r)
    cP$CI.l=tstrsplit(cP$CI.l, '_', type.convert=T)[[2]]
    cP$CI.r=tstrsplit(cP$CI.r, '_', type.convert=T)[[2]]
  
    cP$tchr=as.character(transcript.data$chromosome_name)[match(cP$transcript, transcript.data$gene_id)]
    cP$tpos=transcript.data$start_position[match(cP$transcript, transcript.data$gene_id)]
    cP=cP[cP$tchr %in% rev(chroms), ]#c("X","V","IV","III","II","I"),]
    cP$tchr=factor(cP$tchr,levels=rev(chroms)) #c("X","V","IV","III","II","I"))
    cP$chr=factor(cP$chr,levels=chroms)        #rev(c("X","V","IV","III","II","I")))
    return(cP)
}

plotEQTL=function(cPf, titleme='',CI=T) {
  a=ggplot(cPf,aes(x=pos,y=tpos)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
        geom_point()+
        #,color=cPf$Beta/cPf$SE
        #scale_colour_viridis_c()+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
        facet_grid(tchr~chr,scales="free", space='free')+
        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))+
        ggtitle(titleme)
    if(CI) {
        a=a+geom_segment(aes(x = CI.l, y = tpos, xend = CI.r, yend = tpos)) 
    }
    return(a)
}

plotHotspot=function(cPf, titleme='') {
      h=ggplot(cPf,aes(x=pos)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
        geom_histogram(binwidth=50000)+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
        facet_grid(~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))+
        ggtitle(titleme)

    hd=ggplot_build(h)$data[[1]]
    sig.hp=qpois(1-(.05/length(hd$y)),ceiling(mean(hd$count)))+1 
    h+geom_hline(yintercept=sig.hp)
}

plotHotspot2=function(cPf, titleme='') {
      h=ggplot(cPf,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
        geom_bar(stat="identity", width=50000)+
        xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
        facet_grid(~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))+
        ggtitle(titleme)

       sig.hp=qpois(1-(.05/length(cPf$pos)),ceiling(mean(cPf$count)))+1 
    h+geom_hline(yintercept=sig.hp)
}

makeBinTable=function(cPf, cbin50k) {
    tlinks=makeGRangesFromDataFrame(data.frame(chr=cPf$chr, start=cPf$pos, end=cPf$pos, strand="*")) 

    tlink.cnt=countOverlaps(unlist(cbin50k), tlinks)
    bin.table=data.frame(chr=seqnames(unlist(cbin50k)),
                         pos=round(start(unlist(cbin50k))+25000),
                         count=tlink.cnt)
    return(bin.table)
}

getBinAssignment=function(cPf, cbin50k) {
    tlinks=makeGRangesFromDataFrame(data.frame(chr=cPf$chr, start=cPf$pos, end=cPf$pos, strand="*")) 

    tlink.cnt=findOverlaps(tlinks, unlist(cbin50k))
      
    ub=unlist(cbin50k)
    ub.id=paste0(seqnames(ub), ':', start(ub), '-', start(ub)+width(ub)-1)
    cPf$bin=ub.id[subjectHits(tlink.cnt)]

  
    return(cPf)
}

# given a statistic that increases as a function significance 
# and the same statistic calculated from permutations
# and a table that maps transcripts to the closest marker for each transcript
# calculate functions that map the statistic to FDR
getFDRfx=function(LODr,LODperm, cisMarkers=NULL){

    max.obsLOD=rowMaxs(LODr, value=T) #apply(abs(rff),1,max)
    
    vint=seq(0,max(max.obsLOD)+.01, .01)
    
    ll1=apply(LODperm,3,function(x) rowMaxs(x, value=T))

    # global FDR  ---------------------------------------------------------------
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLOD>thresh) }   )
    names(obsPcnt) = vint

    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(apply(ll1,2,function(x) sum(x>thresh))))
                 })
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    
    pFDR = rev(cummax(rev(pFDR)))
    g.fdrFX=approxfun(pFDR, vint, ties=mean, rule=2)
    g.rtoFDRfx=approxfun(vint,pFDR, ties=mean, rule=2)

    if(!is.null(cisMarkers)) {
        # local FDR --------------------------------------------------------------------
        #ll2=lapply(ll, function(x) x$cmax)
        cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(LODr)),]
        max.obsLODc=LODr[cbind(cMsubset$transcript, cMsubset$marker)]
        ll2=apply(LODperm,3, function(x) x[cbind(cMsubset$transcript, cMsubset$marker)])

        obsPcnt = sapply(vint, function(thresh) { sum(max.obsLODc>thresh) }   )
        names(obsPcnt) = vint
        
        expPcnt = sapply(vint,  
                     function(thresh) { 
                         return(mean(apply(ll2,2,function(x) sum(x>thresh))))
                     })

        names(expPcnt) = vint 
        pFDR = expPcnt/obsPcnt
        
        pFDR[is.na(pFDR)]=0
        pFDR[!is.finite(pFDR)]=0
        #to make sure this is monotonic
        
        pFDR = rev(cummax(rev(pFDR)))
        c.fdrFX=approxfun(pFDR, vint, ties=mean, rule=2)
        c.rtoFDRfx=approxfun(vint,pFDR, ties=mean, rule=2)
        #---------------------------------------------------------------------------------

        return(list(g.fdrFX=g.fdrFX,g.rtoFDRfx=g.rtoFDRfx,
                    c.fdrFX=c.fdrFX, c.rtoFDRfx=c.rtoFDRfx))
    }

    return(list(g.fdrFX=g.fdrFX,g.rtoFDRfx=g.rtoFDRfx))
}




       
chroms=paste0('chr', as.roman(1:16)) 
#c('I','II', 'III', 'IV', 'V', 'X')

#some genome annotation
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))
sgd.genes=getSGD_GeneIntervals(reference.dir)
sgd.genes$chromosome_name=seqnames(sgd.genes)
sgd.genes$start_position=start(sgd.genes)

# for now load all crosses--------------------------------------------------
load(paste0(reference.dir, 'cross.list.RData'))
load(paste0(reference.dir, 'parents.list.RData'))
#reorganize
crossL=list('A'=cross.list[['A']], 
            'B'=cross.list[['B']], 
            '3004'=cross.list[['3004']],
            '393'=cross.list[['393']]

)
parentsL=list('A'=parents.list[['A']], 
              'B'=parents.list[['B']], 
              '3004'=parents.list[['3004']],
              '393'=parents.list[['393']]
)
rm(cross.list)
rm(parents.list)

# need genetic map
genetic.maps=lapply(crossL, extractGeneticMapDf)
#----------------------------------------------------------------------------

# set some variables---------------------------------------------------------
base.dir='/data/single_cell_eQTL/yeast/'
cranger.dir='filtered_feature_bc_matrix/'

#GT.YJM145x : num [1:69385] 0 1 0 1 1 0 0 1 1 2 ...
#  .. ..$ GT.YPS163a

#.$ GT.CBS2888a : num [1:78041] 0 1 1 2 0 0 0 0 0 0 ...
#  .. ..$ GT.YJM981x  : num [1:78041] 1 0 0 0 1 1 1 1 1 1 ...

cList=list(
           '01_2444_44_1-2'='B',
           '02_2444_44-'='B',
           '03_ByxRM_51-_1-2'='A',
           '04_ByxRM_51-'='A',
           '00_BYxRM_480MatA_1'='A',
           '00_BYxRM_480MatA_2'='A',
           '08_2444_cross_10k_Feb_21'='B',
           '09_2444_cross_5k_Feb_21'='B',
           '10_3004_cross_10k_Feb_21'='3004',
           '11_3004_cross_5k_Feb_21'='3004',
           '07_2444-cross-1'='B',
           '07_2444-cross-2'='B',
           '19_393_10k_May10'='393',
           '20_393_20k_May10'='393',
           '21_3004_10k_May10'='3004',
           '22_3004_20k_May10'='3004'
)

hList=list(
           '01_2444_44_1-2'=2^6.1,
           '02_2444_44-'= 2^7.5,
           '03_ByxRM_51-_1-2'=2^7,
           '04_ByxRM_51-'=2^7.4,
           '00_BYxRM_480MatA_1'=2^6.5,
           '00_BYxRM_480MatA_2'=2^6.5,
           '08_2444_cross_10k_Feb_21'=2^5.2,
           '09_2444_cross_5k_Feb_21'= 2^4.9,
           '10_3004_cross_10k_Feb_21'= 2^5.6,
           '11_3004_cross_5k_Feb_21'=  2^4.9,
           '07_2444-cross-1'  = 2^6,
           '07_2444-cross-2'  = 2^5.5,
           '19_393_10k_May10' = 2^5.2,
           '20_393_20k_May10' = 2^5.1,
           '21_3004_10k_May10'= 2^6,
           '22_3004_20k_May10'= 2^6.7
)


experiments=names(cList)
cdesig=as.vector(sapply(cList, function(x) x))
het.thresholds=as.vector(sapply(hList, function(x) x))
data.dirs=paste0(base.dir, 'processed/', experiments, '/')
het.across.cells.threshold=.1

cl <- makeCluster(36)
clusterEvalQ(cl, {
       library(fastglm)
       library(Matrix)
       library(MASS) 
       library(nebula)
       NULL  })

nUMI_thresh=10000
minInformativeCellsPerTranscript=128

nperm=5

#note sgd gene def is updated to include 3' utr               
#sgd.genes=sgd.granges[sgd.granges$type=='gene',]
#sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)

# RUN HMM 

# combine experiments

#good YPSxYJM was 11 and 12
#7,8,

sets=list(
    '3004'=c(9,10),
    'A'=c(3,4),
    'B'=c(1,2,11,12))


sets=list(
          'Ap'=c(5,6)
          )
 
