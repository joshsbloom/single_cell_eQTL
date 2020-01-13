 
generatePlot = function(jointPeaks,dfReturn=F,cutoff=0.1)
{
  colnames(jointPeaks)[colnames(jointPeaks)=="peak.marker"] = "marker"
  jointPeaks = data.frame(jointPeaks,chr=gsub("_+.*","",jointPeaks$marker),pos=as.numeric(gsub(".*_","",jointPeaks$marker)))
  jointPeaksGeneAnnots = geneAnnotsCur[match(jointPeaks$transcript,geneAnnotsCur$wbps_gene_id),]
  jointPeaks = jointPeaks[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  jointPeaksGeneAnnots = jointPeaksGeneAnnots[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  plot_df = data.frame(jointPeaks,jointPeaksGeneAnnots) %>% filter(FDR<cutoff)
  plot_df$chromosome_name = factor(plot_df$chromosome_name,levels=c("X","V","IV","III","II","I"))
  if(dfReturn)
  {
    plot_df
  }else
  {
 (ggplot(plot_df,aes(x=pos,y=start_position)) + facet_grid(chromosome_name~chr,scales="free") + geom_point(aes(alpha=-log10(FDR))) + theme_classic())
  }
}


getFDRfx=function(LODr,LODperm, cisMarkers){

    max.obsLOD=rowMaxs(LODr, value=T) #apply(abs(rff),1,max)
    cMsubset=cisMarkers[which(cisMarkers$transcript %in% rownames(LODr)),]
    
    max.obsLODc=LODr[cbind(cMsubset$transcript, cMsubset$marker)]
    vint=seq(0,max(max.obsLOD)+.05, .05)
    
    ll1=apply(LODperm,3,function(x) rowMaxs(x, value=T))
    ll2=apply(LODperm,3, function(x) x[cbind(cMsubset$transcript, cMsubset$marker)])


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
    pFDR[1]=pFDR[1]-(1e-3)
    g.fdrFX=approxfun(pFDR, vint, ties=mean)
    g.rtoFDRfx=approxfun(vint,pFDR, ties=mean)


    # local FDR --------------------------------------------------------------------
    #ll2=lapply(ll, function(x) x$cmax)
    obsPcnt = sapply(vint, function(thresh) { sum(max.obsLODc>thresh) }   )
    names(obsPcnt) = vint
    
    expPcnt = sapply(vint,  
                 function(thresh) { 
                     return(mean(apply(ll2,2,function(x) sum(x>thresh))))
                 })

    #expPcnt = sapply(vint,  
    #             function(thresh) { 
    #                 return(mean(sapply(ll2,function(x) sum(x>thresh))))
    #             })
    names(expPcnt) = vint 
    pFDR = expPcnt/obsPcnt
    
    pFDR[is.na(pFDR)]=0
    pFDR[!is.finite(pFDR)]=0
    #to make sure this is monotonic
    
    pFDR = rev(cummax(rev(pFDR)))
    pFDR[1]=pFDR[1]-(1e-3)
    c.fdrFX=approxfun(pFDR, vint, ties=mean)
    c.rtoFDRfx=approxfun(vint,pFDR, ties=mean)
    #---------------------------------------------------------------------------------

    return(list(g.fdrFX=g.fdrFX,g.rtoFDRfx=g.rtoFDRfx,
                c.fdrFX=c.fdrFX, c.rtoFDRfx=c.rtoFDRfx))
}



# get peaks from genomewide scan ------------------------------------------------------------
getGlobalPeaks=function(rff, LODr,fdrfx.fc,markerGR,transcript.data, chroms) {
    cPeaksT=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(rff))
        rS=rff[,moi]
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
                             LOD=LODstat, 
                             sbeta=mstat,
                             FDR=fdrfx.fc[['g.rtoFDRfx']](LODstat) )
        cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
        cPeaksT[[cc]]$FDR[(cPeaksT[[cc]]$FDR)>1]=1
    }
 
    cP=rbindlist(cPeaksT, idcol='chrom')
    cP$peak.marker=as.character(cP$peak.marker)
    cP$CI.l=as.character(cP$CI.l)
    cP$CI.r=as.character(cP$CI.r)
    
    cP$peakGpos=markerGR$gcoord[match(cP$peak.marker, colnames(rff))]
    cP$tGpos=transcript.data$gcoord[match(cP$transcript, transcript.data$wormbase_gene)]
    cP$peakLGpos=markerGR$gcoord[match(cP$CI.l, colnames(rff))]
    cP$peakRGpos=markerGR$gcoord[match(cP$CI.r, colnames(rff))]
    #saveRDS(cP, file='/data/single_cell_eQTL/elegans/results/XQTL_F4_2_jointPeaks_filt.RDS')
    #cP=readRDS('/data/single_cell_eQTL/elegans/results/XQTL_F4_2_jointPeaks_filt.RDS')
    cP$tchr=transcript.data$chromosome_name[match(cP$transcript, transcript.data$wormbase_gene)]
    cP$chr=tstrsplit(cP$peak.marker, '_')[[1]]
    cP$pos=as.numeric(tstrsplit(cP$peak.marker, '_')[[2]])
    return(cP)
}
#---------------------------------------------------------------------------------------------------

#tcP_nb=negbinRefit(tcP, Yk, mmp0, fdr.thresh=.1)


negbinNull=function(tnames, Yk, mmp1){
    nbregs=list()  
    for(gn in tnames){
        print(match(gn, colnames(Yk)))
        nbregs[[gn]]=(negbin.reg(Yk[,gn], mmp1[,-1]) )
      }
    return(nbregs)
}
negbinRefitHotspots=function(pms, Yk, Gsub, mmp1, nbM0, fdr.thresh=.1) {
    
    #tcP_nb=tcP[tcP$FDR<fdr.thresh,]
    #tcP_nb$transcript=droplevels(tcP_nb$transcript)
    #tcP_nb=split(tcP_nb, tcP_nb$transcript)

    tcP_nb=list()
    for(gn in colnames(Yk)){
            print(match(gn, colnames(Yk)))
            #nmodelf= negbin.reg(Yk[,gn], mmp1[,-1])     
            nmodelf=nbM0[[gn]]
            #test=fastglm(mmp1[,-1], Yk[,gn], family=poisson(link='log'))
            #test=fastglm(mmp1, Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))

            nmodel=nmodelf$info[3]     
            if(is.na(nmodel)) {next;}
            #pms=as.character(tcP_nb[[gn]]$peak.marker)

            plr=rep(0, length(pms))
            names(plr)=pms
            plb=rep(NA,length(pms))
            names(plb)=pms 
            pls=rep(NA,length(pms))
            names(pls)=pms

            for(p in pms) {
                fm=tryCatch({
                    #negbin.reg(Yk[,gn], cbind(mmp1[,-1],Gsub[,p]))
                    fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))
                }
                    ,error=function(e) {return(NULL);})
                  if(is.null(fm)){
                      plr[p]=NA
                      plb[p]=NA
                      pls[p]=NA
                  } else {
                      msize=length(coef(fm))
                      #plr[p]=fm$info[3]
                      #plb[p]=fm$be[nrow(fm$be),1]
                      plr[p]=as.numeric(logLik(fm)) #fm$info[3]
                      plb[p]=as.numeric(coef(fm)[msize]) #fm$be[nrow(fm$be),1]
                      pls[p]=as.numeric(fm$se[msize]) #fm$be[nrow(fm$be),1]

                  }
            }
            tcP_nb[[gn]]$hotspots=pms
            tcP_nb[[gn]]$negbin.beta=plb
            tcP_nb[[gn]]$negbin.se=pls
            tcP_nb[[gn]]$negbin.LRS=-2*(nmodel-plr)
            tcP_nb[[gn]]$negbin.p=pchisq( tcP_nb[[gn]]$negbin.LRS,1, lower.tail=F) #-2*(nmodel-plr)
            print(gn)
            print(tcP_nb[[gn]])
    }
    return(tcP_nb)
}




negbinRefit=function(tcP, Yk, Gsub, mmp1, nbM0, fdr.thresh=.1) {
    tcP_nb=tcP[tcP$FDR<fdr.thresh,]
    tcP_nb$transcript=droplevels(tcP_nb$transcript)
    tcP_nb=split(tcP_nb, tcP_nb$transcript)


    for(gn in names(tcP_nb)){
            print(match(gn, names(tcP_nb)))
            #nmodelf= negbin.reg(Yk[,gn], mmp1[,-1])     
            nmodelf=nbM0[[gn]]
            #test=fastglm(mmp1[,-1], Yk[,gn], family=poisson(link='log'))
            #test=fastglm(mmp1, Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))

            nmodel=nmodelf$info[3]     
            if(is.na(nmodel)) {next;}
            pms=as.character(tcP_nb[[gn]]$peak.marker)

            plr=rep(0, length(pms))
            names(plr)=pms
            plb=rep(NA,length(pms))
            names(plb)=pms 
            pls=rep(NA,length(pms))
            names(pls)=pms

            for(p in pms) {
                fm=tryCatch({
                    #negbin.reg(Yk[,gn], cbind(mmp1[,-1],Gsub[,p]))
                    fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))
                }
                    ,error=function(e) {return(NULL);})
                  if(is.null(fm)){
                      plr[p]=NA
                      plb[p]=NA
                      pls[p]=NA
                  } else {
                      msize=length(coef(fm))
                      #plr[p]=fm$info[3]
                      #plb[p]=fm$be[nrow(fm$be),1]
                      plr[p]=as.numeric(logLik(fm)) #fm$info[3]
                      plb[p]=as.numeric(coef(fm)[msize]) #fm$be[nrow(fm$be),1]
                      pls[p]=as.numeric(fm$se[msize]) #fm$be[nrow(fm$be),1]

                  }
            }
            tcP_nb[[gn]]$negbin.beta=plb
            tcP_nb[[gn]]$negbin.se=pls
            tcP_nb[[gn]]$negbin.LRS=-2*(nmodel-plr)
            tcP_nb[[gn]]$negbin.p=pchisq( tcP_nb[[gn]]$negbin.LRS,1, lower.tail=F) #-2*(nmodel-plr)
            print(gn)
            print(tcP_nb[[gn]])
    }
    return(tcP_nb)
}

cisNegbinFit=function(cMsubset, Yk, Gsub, mmp1, nbM0) {
    cisNB=list() 
    for(i in 1:nrow(cMsubset)){
         print(i)
         gn=cMsubset[i,'transcript']
         nmodelf=nbM0[[gn]]
         nmodel=nmodelf$info[3]     

         if(is.na(nmodel)) {next;}
         p=cMsubset[i,'marker']
         
         plr=0
         plb=NA
         pls=NA
         fm=tryCatch({
                    #negbin.reg(Yk[,gn], cbind(mmp1[,-1],Gsub[,p]))
                    fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))
                }
                    ,error=function(e) {return(NULL);})
         if(is.null(fm)){
                      plr=0
                      plb=NA
                      pls=NA
         } else {
                      msize=length(coef(fm))
                      #plr[p]=fm$info[3]
                      #plb[p]=fm$be[nrow(fm$be),1]
                      plr=as.numeric(logLik(fm)) #fm$info[3]
                      plb=as.numeric(coef(fm)[msize]) #fm$be[nrow(fm$be),1]
                      pls=as.numeric(fm$se[msize]) #fm$be[nrow(fm$be),1]
                      
         }
        cisNB[[gn]]$negbin.beta=plb
        cisNB[[gn]]$negbin.se=pls
        cisNB[[gn]]$negbin.LRS=-2*(nmodel-plr)
        cisNB[[gn]]$negbin.p=pchisq( cisNB[[gn]]$negbin.LRS,1, lower.tail=F) #-2*(nmodel-plr)
     }
    return(cisNB)

}


getHotspots=function(tcP, gmapd, fdr.thresh=.05, pFDR=T) {

    ####### identify hotspots 
    if(pFDR) {
        cPsig=tcP[tcP$FDR<fdr.thresh,] 
    } else{
        cPsig=tcP[tcP$FDR.negbin<fdr.thresh,]
    }
    
    cPsig.nocis=cPsig[cPsig$tchr!=cPsig$chr,]
    cPsig.nocis$cm5bin=gmapd$cm5bin[match(cPsig.nocis$peak.marker, gmapd$marker)]
    
    bcount=data.frame(values=1:130, lengths=0)

    bcountr=rle(sort(cPsig.nocis$cm5bin))
    bcount$lengths[bcountr$values]=bcountr$lengths
    #bcount=data.frame(values=bcount$values, lengths=bcount$lengths)
    #bcount.holder$lengths[

    gmapdc=split(gmapd, gmapd$chrom)

    sig.hp=qpois(1-(.05/max(gmapd$cm5bin)),ceiling(mean(bcount$lengths)))+1 
    print(sig.hp)
    gcm=t(sapply(gmapdc, function(x) range(x$cm5bin)))
    bcount$chr=rep(chroms, (1+gcm[,2]-gcm[,1]))
    bcount$sig=bcount$lengths>sig.hp

    bcounts=split(bcount,bcount$chr)

    #automate peak detection
    bcl=lapply(bcounts, function(x) {
           f=findpeaks(x$lengths,minpeakheight=sig.hp, minpeakdistance=3, nups=0)[,2]
           x$peaks=FALSE
           x$peaks[f]=T
           x$peaks[x$sig==F]=F
           return(x)  })
    hotspot.bins=unlist(sapply(bcl, function(x) x$values[x$peaks]))
    cpb=cPsig.nocis[cPsig.nocis$cm5bin %in% as.vector(unlist(hotspot.bins)), ]
    cpbs=split(cpb, cpb$cm5bin)
    hotspot.markers=as.character(sapply(cpbs, function(x) names(which.max(table(x$peak.marker)))))
    return(list(hotspot.markers=hotspot.markers, bin.counts=bcl))
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

# bin markers to genetic map 
# bin. size is centimorgans
binMarkersToGeneticMap=function(Gr, gmap, bin.size=5) {
    gmapd=gmap[match(colnames(Gr), gmap$marker),]
    gmapdc=split(gmapd, gmapd$chrom)
    gmap.max=sapply(gmapdc, function(x) ceiling(max(x$map)))
    mbin=0
    #bin.size=5
    for(cc in chroms) {
        fib=findInterval(gmapdc[[cc]]$map,seq(0,gmap.max[[cc]]+(bin.size-gmap.max[[cc]]%%bin.size),bin.size))
        gmapdc[[cc]]$cm5bin=fib+mbin
        mbin=max( gmapdc[[cc]]$cm5bin )
    }
    gmapd=rbindlist(gmapdc)
    return(gmapd)
}


    # no point keeping all this, just retain max stats 
    #ll=replicate(nperm, {
    #                 nind=sample(1:nrow(Yfe))
    #                pme=abs(crossprod(Yfe[nind,],Gf)/(nrow(Yfe)-1))
    #                lout=list(amax=rowMaxs(pme, value=T),
    #                          cmax=pme[cbind(cMsubset$transcript, cMsubset$marker)] 
    #                )
    #                return(lout)
    #    }, simplify=F )

    #fix this ---- 
#    llM=-nrow(Yfe)*log(1-ll$amax^2)/(2*log(10)) 
#    llC=-nrow(Yfe)*log(1-ll$cmax^2)/(2*log(10))
# ------------
