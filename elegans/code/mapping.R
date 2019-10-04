getFDRfx=function(r, Yre,Gr,nperm=5, vint=seq(0.001,.15,.001)){
    
    max.obsLOD=apply(abs(r),1,max)
    # no point keeping all this, just retain max stats 
    ll=replicate(nperm, {
                    nind=sample(1:nrow(Yre))
                    rowMaxs(abs(crossprod(Yre[nind,],Gr)/(nrow(Yre)-1)),value=T)
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
    fdrFX=approxfun(pFDR, vint, ties=mean)
    rtoFDRfx=approxfun(vint,pFDR, ties=mean)

    return(list(fdrFX=fdrFX,rtoFDRfx=rtoFDRfx))
}


# cells with NAs are useless
crit=!is.na(fine.classification$broad) & !is.na(fine.classification$fine) 
# & log2(fine.classification$total)>8.5 & log2(fine.classification$total)<15.5

joint.covariates=data.frame(fine.classification)
joint.covariates=joint.covariates[crit,]

Y=t(as.matrix(counts))
Yr=Y[crit,]

#calc Variance for each transcript, if zero, boot it
Yr.var=colVars(Yr, parallel=T)
Yr=Yr[,Yr.var>0]

Gr=geno[crit,]

G=standardise2(Gr)
colnames(G)=colnames(G)
rownames(G)=rownames(Yr)

#model total reads, batch, broad and fine tissue classifications 
mm=model.matrix(Yr[,1]~log2(joint.covariates$total)+joint.covariates$Batch+joint.covariates$broad_tissue) #+joint.covariates$fine_tissue)
BOLS=lm.fit(mm, log2(Yr+1), intercept=F)

Yresid=residuals(BOLS)
rownames(Yresid)=rownames(Yr)
colnames(Yresid)=colnames(Yr)
rm(BOLS)
#Yresid.var=colVars(Yresid,parallel=T)

Yre=standardise2(Yresid)
rownames(Yre)=rownames(Yresid)
colnames(Yre)=colnames(Yresid)

# identify closest marker to transcript for cis only test (to do) --------------------------------------------
g2m=tstrsplit(colnames(G),'_')
markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                             start=as.numeric(g2m[[2]]),
                                             end=as.numeric(g2m[[2]]), 
                                             strand="*",name=colnames(G)))


tag=transcript.annotations[transcript.annotations$type=='gene',]
tag=tag[match(colnames(Yre), tag$gene_id),]
cisMarkers=data.frame(transcript=colnames(Yre), marker=colnames(G)[nearest(tag,markerGR)],stringsAsFactors=F)
cisMarkers=na.omit(cisMarkers)
#indices of cis eQTL in 
#rows are transcripts, columns are markers
cis.index=cbind(match(cisMarkers$transcript, colnames(Yre)), match(cisMarkers$marker, colnames(G)))
#------------------------------------------------------------------------------------------------------------


# 1D scan speedup 
r=crossprod(Yre,G)/(nrow(Yre)-1)
LOD=-nrow(Yre)*log(1-r^2)/(2*log(10))

#pool for permutations
wgFDR= getFDRfx(r, Yre, G, nperm=5)

#get peaks and confidence intervals------------------------------
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
#---------------------------------------------------------------
cP=do.call('rbind', cPeaks)
plot(match(cP$peak.marker[cP$FDR<.01], colnames(G)),
     match(cP$transcript[cP$FDR<.01], transcript.data$wormbase_gene), 
      xlab='marker index', ylab='transcript index', main='joint analysis FDR < 1%')


# Map within each tissue type --------------------------------------------------
qtl.maps=list()
qtl.fdrs=list()
rho.matrices=list()
residual.matrices=list()
fc=unique(as.character(joint.covariates$fine_tissue))
fc=fc[!is.na(fc)]
bc=unique(as.character(joint.covariates$broad_tissue))
bc=bc[!is.na(bc)]
types=c(bc,fc)
cty=c(rep('broad', length(bc)), rep('fine', length(fc)))

for(kk in 1:length(types)) {
    uc=types[kk]
    print(uc)
    if(cty[kk]=='broad') {
        tk=which(joint.covariates$broad_tissue==uc)
        print('broad')
    }
    if(cty[kk]=='fine') {
        tk=which(joint.covariates$fine_tissue==uc)
        print('fine')
    }
    
    covs=joint.covariates[tk,]
    Yk=Yr[tk,]
    Gsub=Gr[tk,]

    # OLS, replace with glm for next iteration
    # presumably run once only
    BOLS=lm(log2(Yk+1)~log2(covs$total)+covs$Batch)
    Yresid=residuals(BOLS)
    rownames(Yresid)=rownames(Yk)
    colnames(Yresid)=colnames(Yk)
    rm(BOLS)
    residual.matrices[[uc]]=Yresid

    #so slow ...
    Yfe=standardise2(Yresid)
    rownames(Yfe)=rownames(Yresid)
    colnames(Yfe)=colnames(Yresid)
    Gf=standardise2(Gsub)
    colnames(Gf)=colnames(G)
    rownames(Gf)=rownames(Yfe)
    #transcripts that are singular after the regression
    shitT=unique(which(is.na(Yfe), arr.ind=T)[,2])
    Yfe=Yfe[,-shitT]

    # then map within each classification
    rff=crossprod(Yfe,Gf)/(nrow(Yfe)-1)
    rho.matrices[[uc]]=rff
    LODr=-nrow(Yfe)*log(1-rff^2)/(2*log(10))

    fdrfx.fc=getFDRfx(rff,Yfe,Gf,nperm=5,vint=seq(0.001,.45,.001))
    qtl.fdrs[[uc]]=fdrfx.fc
   
    cPeaksT=list()
    for(cc in chroms) {
        print(cc)
        moi=grep(paste0('^',cc,'_'), colnames(rff))
        rS=rff[,moi]
        #mstat=apply(rS,1,max)
        mstat.pos=apply(abs(rS),1,which.max)
        lookup=cbind(1:length(mstat.pos), mstat.pos)
        mstat=rS[lookup]
        LODstat=LODr[,moi][lookup]
        CIs=matrix(NA,length(mstat),2)
        cmoi=colnames(LODr[,moi])
        for(peak in 1:length(mstat)){
            CIs[peak,]=cmoi[range(which(LODr[peak,moi]>LODstat[peak]-1.5))]
        }   

        cPeaksT[[cc]]=data.frame( transcript=rownames(rS),
                             peak.marker=colnames(rS)[mstat.pos],
                             CI.l=CIs[,1],
                             CI.r=CIs[,2],
                             LOD=LODstat, 
                             sbeta=mstat,
                             FDR=fdrfx.fc[[2]](abs(mstat)) )
        cPeaksT[[cc]]$FDR[is.na(cPeaksT[[cc]]$FDR)]=1
    }
 
   qtl.maps[[uc]]=rbindlist(cPeaksT, idcol='chrom')
   print(sum(qtl.maps[[uc]]$FDR<.2, na.rm=T))
}
#---------------------------------------------------------------------------------------------------






cPt=do.call('rbind',cPeaksT)

plot(match(cPt$peak.marker[cPt$FDR<.9], colnames(G)),
     match(cPt$transcript[cPt$FDR<.9], transcript.data$wormbase_gene), 
      xlab='marker index', ylab='transcript index', main='joint analysis FDR < 1%')
















































#new genotypes 
r2=crossprod(Yre,Gr1b)/(nrow(Yre)-1)


cis1=-nrow(Yre)*log(1-r1[cis.index]^2)/(2*log(10))
cis2=-nrow(Yre)*log(1-r2[cis.index]^2)/(2*log(10))



cisP_G=rep(NA, nrow(cisMarkers))
cisP_G2=rep(NA, nrow(cisMarkers))
for(rr in 1:nrow(cisMarkers)) {
    print(rr)
    lmm=lm(Yre[,cisMarkers$transcript[rr]]~G2[crit,cisMarkers$marker[rr]])
    plot(G2[crit,cisMarkers$marker[rr]], Yre[,cisMarkers$transcript[rr]],col="#00000022") 
    abline(lm(Yre[,cisMarkers$transcript[rr]]~G2[crit,cisMarkers$marker[rr]]), col='red')
    plot(G2[crit,cisMarkers$marker[rr]], Y[crit,cisMarkers$transcript[rr]],col="#00000022") 

    cisMarker=G2[crit,cisMarkers$marker[rr]]
    cisTranscriptL=log2(Y[crit,cisMarkers$transcript[rr]]+1)
    cisTranscriptR=(Y[crit,cisMarkers$transcript[rr]])

    cisTissues=joint.covariates$broad_tissue
    lmm=lm(cisTranscriptL~log2(joint.covariates$total)+joint.covariates$Batch+cisMarker*cisTissues)

    glmm=glm.nb(cisTranscriptR~offset(log(joint.covariates$total))+joint.covariates$Batch+cisMarker*cisTissues)
    
    
    G2[crit,cisMarkers$marker[rr]])


    cisP_G2[rr]=(drop1(lmm, test='Chisq'))[[5]][2]
    lmm=lm(Yre[,cisMarkers$transcript[rr]]~rank(G[crit,cisMarkers$marker[rr]]) )
    cisP_G[rr]=(drop1(lmm, test='Chisq'))[[5]][2]
}












# why is this so fucking slow 
Gr1=standardise2(colRanks(G[crit,],parallel=T))
colnames(Gr1)=colnames(G)
rownames(Gr1)=rownames(Yre)

r=crossprod(Yre,Gr1)/(nrow(Yre)-1)
LOD1=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL1=which(LOD1>8, arr.ind=T)
rm(r)
rm(LOD1)
plot(sigL1[,2], sigL1[,1], pch=21,main='ranks r/qtl genotypes')
#------------------

#-------------------------------------------------
Gr2=standardise2(colRanks(G2[crit,],parallel=T))
colnames(Gr2)=colnames(G)
rownames(Gr2)=rownames(Yre)

r=crossprod(Yre,Gr2)/(nrow(Yre)-1)
LOD2=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL2=which(LOD2>8, arr.ind=T)
plot(sigL2[,2], sigL2[,1], pch=21,main='ranks custom genotypes')

rm(r)
rm(LOD2)
#-------------------------------------------------


#-------------------------------------------------
Gr2b=standardise2(G2[crit,])
colnames(Gr2b)=colnames(G)
rownames(Gr2b)=rownames(Yre)
r=crossprod(Yre,Gr2b)/(nrow(Yre)-1)
LOD3=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL3=which(LOD3>8, arr.ind=T)
plot(sigL3[,2], sigL3[,1], pch=21,main='custom genotypes')

rm(r)
rm(LOD3)
#---------------------------------------------------

#------------------------------------------
Gr1b=standardise2(G[crit,])
colnames(Gr1b)=colnames(G)
rownames(Gr1b)=rownames(Yre)
r=crossprod(Yre,Gr1b)/(nrow(Yre)-1)
LOD4=-nrow(Yre)*log(1-r^2)/(2*log(10))
sigL4=which(LOD4>8, arr.ind=T)
rm(r)
rm(LOD4)
x11()
plot(sigL4[,2], sigL4[,1], pch=21,main='r/qtl genotypes')


par(mfrow=c(2,2))
plot(sigL1[,2], sigL1[,1], pch=21, cex=.1, col="#00000022", main='ranks r/qtl genotypes',xlab='marker index', ylab='transcript index')
plot(sigL4[,2], sigL4[,1], pch=21,, cex=.1, col="#00000022", main='r/qtl genotypes',xlab='marker index', ylab='transcript index')
plot(sigL2[,2], sigL2[,1], pch=21,, cex=.1, col="#00000022", main='ranks custom genotypes',xlab='marker index', ylab='transcript index')
plot(sigL3[,2], sigL3[,1], pch=21,, cex=.1, col="#00000022", main='custom genotypes',xlab='marker index', ylab='transcript index')



G2sum=colSums(G2)
G2rs=rowSums(G2)



tt=eMats/colSums(eMats)
test=eLikeHMM(tt,ttMat)


#CB=N2
#n=5
#k=1
#nmk=n-k
#bco=choose(n,k)
#ee=.0001

#bco*((1-ee)^k)*(ee^(nmk))
#bco*(.5)^n
#bco*((1-ee)^nmk)*ee^k

#ref/het/alt
#pRR  = dbinom(n-k,n,ee) *      (1-r)/2
#pHet = choose(n,k)*(1/(2^n))*   r
#pAA  = dbinom(k,n,ee) *        (1-r)/2


nind=nrow(cross2$cross_info)
cc=cut(1:nind, 100)

f='1'
gps=readRDS(paste0('/data/single_cell_eQTL/elegans/XQTL_F4_2/hmm_v2/', f))

#501 is weird

for(h in 600:691){
    i=rownames(gps[[1]])[h]# 691
    coii='II'
    coi=paste0('^', coii, '_') #'^II_'

    #plot(gps[[2]][i,1,])
    N_bs=N2.counts[grep(coi, rownames(N2.counts)),i]
    C_bs=CB.counts[grep(coi, rownames(CB.counts)),i]

   # N_bs0=N2.counts0[grep(coi, rownames(N2.counts0)),i]
   # C_bs0=CB.counts0[grep(coi, rownames(CB.counts0)),i]

   # N_bsR=N2.countsR[grep(coi, rownames(N2.countsR)),i]
   # C_bsR=CB.countsR[grep(coi, rownames(CB.countsR)),i]

    cGcall=cross2$geno[[coii]][i,]

    par(mfrow=c(2,1))
    plot(.5*gps[[coii]][i,2,]+gps[[coii]][i,3,],main='p(C)', col='blue', type='l')
    #points(gps[[coii]][i,2,],main='p(NC)', col='grey')
    #points(gps[[coii]][i,3,],main='p(CC)', col='blue')
    #test=unlist(cGcall, function(x) switch(x, 'blue', 'grey', 'orange', 'purple', 'green'))
    cme=which(N_bs>0)
    if(length(cme)==0) {next;}

    plot(cme,N_bs[cme], type='h', ylim=c(-10,10), main='br, strict',col=cGcall[names(cme)]) #sapply(cGcall, function(x) switch(x, 'blue', 'grey', 'orange', 'purple', 'green')))
    points(cme, -C_bs[cme], type='h', col=cGcall[names(cme)]) #, col=cGcall)
    abline(h=0)

   # plot(-N_bs0, type='h', ylim=c(-10,10), main='br')
   # points(C_bs0, type='h') #, col=cGcall)

    #plot(-N_bsR, type='h', ylim=c(-10,10), main='raw')
    #points(C_bsR, type='h') #, col=cGcall)

    print(h)
    readline()
    }
#black is CC
#red is NC
#green is NN
#blue is NN or NC
#cyan is NC or CC





hets=which(cGcall==2)
points(hets, rep(0,length(hets)), col='red', pch=17, cex=2)
notNN=which(cGcall==5)
points(notNN, rep(0,length(notNN)), col='cyan',pch=17,cex=2)
notCC=which(cGcall==4)
points(notCC, rep(0,length(notCC)), col='blue', pch=17, cex=2)
isCC=which(cGcall==1)
points(isCC, rep(0,length(isCC)), col='black', pch=17, cex=2)
isNN=which(cGcall==3)
points(isNN, rep(0,length(isNN)), col='green',pch=17, cex=2)













