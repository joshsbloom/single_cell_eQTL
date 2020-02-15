#elegans chromosomes 
chroms=c('I', 'II', 'III', 'IV', 'V', 'X')

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

# wrapper for abind commmands
acomb1 <- function(...) abind(..., along=1)


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



# get peaks from genomewide scan ------------------------------------------------------------
# also restructure table for ease of plotting 
getGlobalPeaks=function( LODr,fdrfx.fc,markerGR,transcript.data, chroms) {
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
                             LOD=LODstat, 
                             #sbeta=mstat,
                             FDR=fdrfx.fc[['g.rtoFDRfx']](LODstat) )
        cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
        cPeaksT[[cc]]$FDR[(cPeaksT[[cc]]$FDR)>1]=1
    }
 
    cP=rbindlist(cPeaksT, idcol='chrom')
    cP$peak.marker=as.character(cP$peak.marker)
    cP$chr=tstrsplit(cP$peak.marker, '_')[[1]]
    cP$pos=as.numeric(tstrsplit(cP$peak.marker, '_')[[2]])
    cP$CI.l=as.character(cP$CI.l)
    cP$CI.r=as.character(cP$CI.r)
    cP$CI.l=tstrsplit(cP$CI.l, '_', type.convert=T)[[2]]
    cP$CI.r=tstrsplit(cP$CI.r, '_', type.convert=T)[[2]]
  
    cP$tchr=transcript.data$chromosome_name[match(cP$transcript, transcript.data$wormbase_gene)]
    cP$tpos=transcript.data$start_position[match(cP$transcript, transcript.data$wormbase_gene)]
    cP=cP[cP$tchr %in% c("X","V","IV","III","II","I"),]
    cP$tchr=factor(cP$tchr,levels=c("X","V","IV","III","II","I"))
    cP$chr=factor(cP$chr,levels=rev(c("X","V","IV","III","II","I")))
    return(cP)
}
#---------------------------------------------------------------------------------------------------



# formula for calculating the negative binomial glm log likelihood 
nbLL=function (y, mu, theta)  {
    return( -sum ( (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta + y) ))
}

# fast calculation of negative binomial LL across the genome if theta is known 
domap=function(gn, ...) { 
       theta.est=thetas[gn]
       YY=Y[,gn]
       nbn=negative.binomial(theta=theta.est,link='log')
      
       fnbrN=fastglmPure(DM, YY, family=nbn)
       nmLLik=nbLL(fnbrN$y,fnbrN$fitted.value, theta.est)

       LRS=rep(NA,ncol(Gsub))
       names(LRS)=colnames(Gsub)

       XXn=cbind(DM,1)
       idx=1:ncol(Gsub)
       for(gx in idx) { 
           #colnames(Gsub)){
           XXn[,7]=Gsub[,gx]
           fnbrF=fastglmPure(XXn, YY,  family=nbn)
           LRS[gx]=-2*(nmLLik-nbLL(fnbrF$y,fnbrF$fitted.value, theta.est))
        }
       return(LRS)
}


# fit a negative binomial model for each transcript including covariates and the closest 
# marker to each transcript
# use Rfast2::negbin.reg for computational efficiency
# for very large theta estimates or model convergence failures
# use mgcv::gam for more robust algorithm
#gsubs=scale(Gsub)
#sg=svd(gsubs)
#um=sg$u
#colnames(um)=paste0('X',seq(1,ncol(um)))

#try not adding a size parameter 
#test=glmpca(t(Yk), L=10, fam='nb', nb_theta=.1, X=m0)
#emus=sapply(cisNB, function(x) x$emu)
#disps=sapply(cisNB, function(x) x$theta)
#ps=sapply(cisNB, function(x) x$negbin.p)
#
#dffe=data.frame(emus=emus, disps=disps, p=ps)
#dffe=dffe[-which(dffe$disps>100),]

fitNBCisModel=function(cMsubset, mmp1, Yk, Gsub) {
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
        #print(cisNB[[gn]])
    }

   # nbR.var=colVars(nbR)
   # nbR=nbR[,-which(is.na(nbR.var))]
    return(cisNB) #, nbR=nbR))
}

# experimental version using fixest::fenegbin()
fitNBCisModelFE=function(cMsubset, mmp1, Yk, Gsub) {
    mmp0=mmp1[,-1]

    #nbR=Yk
    #nbR[is.numeric(nbR)]=NA
    ##nbP=nbR

    dff=data.frame(batch=covs$Batch,ltotal=log(covs$total))
    m0=model.matrix(dff$ltotal~dff$batch)
    nmsize=ncol(m0)+1
    cisNB=list() 
    for(r in 1:nrow(cMsubset)){
        print(r)
        gn=as.character(cMsubset$transcript[r])
        p=as.character(cMsubset$marker[r])
        #cmodel=cbind(mmp0, Gsub[,p])
        #nbr=negbin.reg(Yk[,gn], cmodel,maxiters=500)
        dff$cis=Gsub[,p]
        dff$Y=Yk[,gn]
        #fnbrF=feglm(Y~batch+cis, offset=dff$ltotal,data=dff,family='poisson',nthreads=24)
        #fnbrN=feglm(Y~batch, offset=dff$ltotal,data=dff,family='poisson',nthreads=24)
        
        fnbrF=fenegbin(Y~batch+cis, offset=dff$ltotal,data=dff,nthreads=24)
        # dffT=cbind(dff, um)
        # fnbrT= update(fnbrF, ~.+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25, data=dffT)
        coef.est=coef(fnbrF)[1:nmsize]
        theta.est=coef(fnbrF)[nmsize+1]
        # XX=cbind(m0, Gsub[,p])
        # fnbrF=fastglm(XX, Yk[,gn], start=coef.est, offset=dff$ltotal, family=negative.binomial(theta=theta.est,link='log'), maxit=500)
        fnbrN=fastglm(m0,Yk[,gn], start=coef.est[-nmsize],  offset=dff$ltotal, family=negative.binomial(theta=theta.est,link='log'), maxit=100)
        cisNB[[gn]]$transcript=gn
        cisNB[[gn]]$lmarker=p

        cisNB[[gn]]$theta=theta.est
        cisNB[[gn]]$emu=mean(predict(fnbrF))

        cisNB[[gn]]$negbin.beta=as.numeric(coef.est[nmsize])
        cisNB[[gn]]$negbin.se=as.numeric(fnbrF$se[nmsize])
        cisNB[[gn]]$LLik=as.numeric(logLik(fnbrF))
        cisNB[[gn]]$negbin.LRS=as.numeric(-2*(logLik(fnbrN)-logLik(fnbrF)))
        cisNB[[gn]]$negbin.p=pchisq( cisNB[[gn]]$negbin.LRS,1, lower.tail=F) 
        #-2*(nmodel-plr)
        cisNB[[gn]]$fmodelBs=coef(fnbrF)
    }
    return(cisNB) 
}



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

# calculate hotspots given a binning of genetic map to 
getHotspots=function(tcP, gmapd, fdr.thresh=.05) {

    ####### identify hotspots 
    cPsig=tcP[tcP$FDR<fdr.thresh,] 
    cPsig.nocis=cPsig[cPsig$tchr!=cPsig$chr,]
    cPsig.nocis$bin=gmapd$cm5bin[match(cPsig.nocis$peak.marker, gmapd$marker)]
    
    bcount=data.frame(values=1:max(gmapd$cm5bin), lengths=0)

    bcountr=rle(sort(cPsig.nocis$bin))
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
    cpb=cPsig.nocis[cPsig.nocis$bin %in% as.vector(unlist(hotspot.bins)), ]
    cpbs=split(cpb, cpb$bin)
    hotspot.markers=as.character(sapply(cpbs, function(x) names(which.max(table(x$peak.marker)))))
    return(list(hotspot.markers=hotspot.markers,bin.counts=bcl, hotspot.linkages=cpbs))
}

# note use performance::r2() for pseudo r^2, convert fastglm to glm  



getMeff_Li_and_Ji=function(cor.mat) {
    evals = eigen(cor.mat,symmetric=T)$values
    M = length(evals)
    L = M-1
    # Equation 5 from Li 
    intevals=ifelse(evals>=1, 1, 0)
    # modification for negative eigenvalues JB
    # add up the contribution of fractional eigenvalues 
    nonintevals=evals[evals>0]-floor(evals[evals>0]) #[evals>0]
    Meff.li=sum(intevals)+sum(nonintevals)
    print(Meff.li)
    return(Meff.li)
}

getFDRmeff=function(pvec,q=.05,Meff) {
#find  the max i, where observed p-value is less than 
pvec < (q/Meff+(((i-1)/(M-1))*(q-(q/Meff))))
}

#functions for computing multivariate scan 
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












# THESE FUNCTIONS BELOW ARE NO LONGER BEING USED ....========================================================


#tcP_nb=negbinRefit(tcP, Yk, mmp0, fdr.thresh=.1)
negbinNullDev=function(tnames, Yk, mmp1, nbM0, covs, threshold=100) { #0,m.to.theta ) {
    #, smoothed.disp){
    logLikV=rep(NA, length(tnames))
    names(logLikV)=tnames
    thetaV=rep(NA, length(tnames))
    names(thetaV)=tnames
    nbR=Yk
    nbR[is.numeric(nbR)]=NA
    mexp=colMeans(Yk)
    for(gn in tnames){
        print(match(gn, tnames))
        nmodelf=nbM0[[gn]]
        theta=nmodelf$info[4]
        if(is.na(theta) | theta>threshold ){
              test=mgcv::gam(Yk[,gn]~covs$Batch, offset=log(covs$total), family=mgcv::nb , method="ML")
              theta=test$family$getTheta(T)
        }
        #theta.smoothed=m.to.theta(mexp[gn])
        #theta.smoothed=smoothed.disp[gn]
        #nbr=mgcv::gam(Yk[,gn]~covs$Batch, offset=log(covs$total), family=mgcv::nb , method="ML")
        nbr= tryCatch({
            #fastglm(mmp1, Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))},
            fastglm(mmp1, Yk[,gn], family=negative.binomial(theta=theta,link='log'))},
            error=function(e){return(NULL);})
        if(is.null(nbr)) {next;}
        logLikV[gn]=logLik(nbr)
        #thetaV[gn]=nmodelf$info[4] #nbr$family$getTheta(T)
        thetaV[gn]=theta
        nbR[,gn]=residuals(nbr,'deviance') 
    }
    Yfe.var=colVars(nbR)
    nbR=nbR[,which(!is.na(Yfe.var))]
    return(list(nbR=nbR, logLik=logLikV, theta=thetaV))
    #return(nbR)
}


negbinNull=function(tnames, Yk, mmp1){
    nbregs=list()  
    for(gn in tnames){
        print(match(gn, colnames(Yk)))
        nbregs[[gn]]=(negbin.reg(Yk[,gn], mmp1[,-1]) )
      }
    return(nbregs)
}

#use theta.ml trick and poisson regression for speedup here
negbinNullP=function(tnames, Yk, mmp2){
    nbregs=list()  
    for(gn in tnames){
        print(match(gn, colnames(Yk)))
       # test=mgcv::gam(Yk[,gn]~mmp1[,-1], family=mgcv::nb , method="ML")

      #  test2=fastglm(mmp1,Yk[,gn],family=poisson(link='log'), method=3)
      #  pt2=predict(test2,mmp1)
      #  nbregs[[gn]]=theta.ml(Yk[,gn], exp(pt2), limit=50)
        nbregs[[gn]]=(negbin.reg(Yk[,gn], mmp2[,-1],maxiters=500) )
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




negbinRefit=function(tcP, Yk, Gsub, mmp1, nbM0, fdr.thresh=.1,sfit) {
    tcP_nb=tcP[tcP$FDR<fdr.thresh,]
    tcP_nb$transcript=droplevels(tcP_nb$transcript)
    tcP_nb=split(tcP_nb, tcP_nb$transcript)


    for(gn in names(tcP_nb)){
            print(match(gn, names(tcP_nb)))
            #nmodelf= negbin.reg(Yk[,gn], mmp1[,-1])     
            nmodelf=nbM0[[gn]]
            #test=fastglm(mmp1[,-1], Yk[,gn], family=poisson(link='log'))
            #test=fastglm(mmp1, Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))

            #nmodelLL=nmodelf$info[3]     
            nmodelLL=sfit$logLik[gn]
            #theta.est=nmodelf$info[4]
            theta.est=sfit$theta[gn]
            
            if(is.na(nmodelLL)) {next;}
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
                    fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=theta.est,link='log'))
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
            tcP_nb[[gn]]$negbin.LRS=-2*(nmodelLL-plr)
            tcP_nb[[gn]]$negbin.p=pchisq( tcP_nb[[gn]]$negbin.LRS,1, lower.tail=F) #-2*(nmodel-plr)
            print(gn)
            print(tcP_nb[[gn]])
    }
    return(tcP_nb)
}

cisNegbinFit=function(cMsubset, Yk, Gsub, mmp1, nbM0, sfit) {
    cisNB=list() 
    for(i in 1:nrow(cMsubset)){
         print(i)
         gn=cMsubset[i,'transcript']
         nmodelf=nbM0[[gn]]
         nmodel=nmodelf$info[3]     
          
         nmodelLL=sfit$logLik[gn]
         #theta.est=nmodelf$info[4]
         theta.est=sfit$theta[gn]

         if(is.na(nmodel)) {next;}
         p=cMsubset[i,'marker']
         
         plr=0
         plb=NA
         pls=NA
         fm=tryCatch({
                    #negbin.reg(Yk[,gn], cbind(mmp1[,-1],Gsub[,p]))
                    #fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=nmodelf$info[4],link='log'))
                    fastglm(cbind(mmp1,Gsub[,p]), Yk[,gn], family=negative.binomial(theta=theta.est,link='log'))

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
        cisNB[[gn]]$negbin.LRS=-2*(nmodelLL-plr)
        cisNB[[gn]]$negbin.p=pchisq( cisNB[[gn]]$negbin.LRS,1, lower.tail=F) #-2*(nmodel-plr)
     }
    return(cisNB)

}



