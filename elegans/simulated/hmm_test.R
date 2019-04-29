library(qtl)
library(qtl2)
library(qtl2convert)
library(data.table)
library(Rfast)
library(abind)
source('/data/single_cell_eQTL/elegans/elegans_eQTL_fx.R')


snps=readRDS('scRNAseqSNPlist.Apr2019.rds')
snp.df=data.frame(markerName=snps$markerNames, chr=snps$markerChrom, 
                  pos=snps$markerPos,
                  map=snps$map, stringsAsFactors=F)

snpsU=readRDS('scRNAseqSNPlist.uniform.Apr2019.rds')
snpU.df=data.frame(markerName=snpsU$markerNames, chr=snpsU$markerChrom, 
                  pos=snpsU$markerPos,
                  map=snpsU$map, stringsAsFactors=F)

wt.power=list()
for(j in 6:10) {
    print(j)
    sim.data=readRDS(paste0('/data/single_cell_eQTL/elegans/simulated/','scRNAseqSimulation_generation_', j,  '_wtMap.rds'))
    real=sim.data[[3]]
    sparse=sim.data[[4]]

    colnames(sparse)=as.character(1:ncol(sparse))
    rownames(sparse)=snp.df$markerName

    gm=snp.df$map
    names(gm)=snp.df$markerName
    gmap.s=split(gm, snp.df$chr)
    cross=buildCrossObject(sparse,gmap.s,'f2','_')
    gps=calc_genoprob(cross, error_prob=.001)

    gX=do.call('abind', gps)
    sim.sites=seq(1,dim(gX)[3],10)
    ps=rep(NA, length(sim.sites))
    for(i in 1:length(sim.sites)){
        #print(i)
        x=real[sim.sites[i],]*2
        y=x*.4+rnorm(1000,.6)
        #lm(y~x)
        ps[i]=anova(lm(y~gX[,,sim.sites[i]]), lm(y~1))[[6]][2]
    }
    wt.power[[j]]=sum(ps<.05/1e3)/length(sim.sites)
}

#visualize

x11()
for( i in 1:100) {
    ind=i

    par(mfrow=c(2,1))
    plot(-1+real[,ind]*2, type='h', ylab='simulated genotype', main= sum(!is.na(sparse[,ind])))
    points(-1+2*sparse[,ind], pch=20,cex=2, col='yellow')
    #points((2*gpi[3,]),col='blue', ylab='p(0)')
    #points((2*(-1*gpi[1,])),col='red', ylab='p(2)')
    #points(2*gpi[2,]-1,col='grey', ylab='p(1)')
    abline(v=cumsum(sapply(gmap.s, length)),lty=2)

    plot(gX[i,3,],col='blue', ylim=c(0,1), ylab='genotype probability')
    points(gX[i,2,],col='grey')
    points(gX[i,1,],col='red')
    points(sparse[,ind],cex=2, pch=20, col='yellow' )
    abline(v=cumsum(sapply(gmap.s, length)),lty=2)

    readline()
}





rec1.power=list()
for(j in 6:10) {
    print(j)
    sim.data=readRDS(paste0('/data/single_cell_eQTL/elegans/simulated/','scRNAseqSimulation_generation_', j,  '_uniformMap.rds'))
    real=sim.data[[3]]
    sparse=sim.data[[4]]

    colnames(sparse)=as.character(1:ncol(sparse))
    rownames(sparse)=snp.df$markerName

    gm=snpU.df$map
    names(gm)=snp.df$markerName
    gmap.s=split(gm, snpU.df$chr)
    
    cross=buildCrossObject(sparse,gmap.s,'f2','_')
    gps=calc_genoprob(cross, error_prob=.001)

    gX=do.call('abind', gps)
    sim.sites=seq(1,dim(gX)[3],10)
    ps=rep(NA, length(sim.sites))
    for(i in 1:length(sim.sites)){
        #print(i)
        x=real[sim.sites[i],]*2
        y=x*.4+rnorm(1000,.6)
        #lm(y~x)
        ps[i]=anova(lm(y~gX[,,sim.sites[i]]), lm(y~1))[[6]][2]
    }
    rec1.power[[j]]=sum(ps<.05/1e3)/length(sim.sites)
}







real=readRDS('miSeqSimulationFull.rds')
r0=apply(real, 1, function(x) sum(x==0))
r1=apply(real, 1, function(x) sum(x==0.5))
r2=apply(real, 1, function(x) sum(x==1))
r.af=(2*r0+r1)/(2*r0+2*r1+2*r2)

sparse=readRDS('miSeqSimulationSparse.rds')
graw=sparse[,1:100]
rownames(graw)=snp.df$markerName
colnames(graw)=as.character(1:ncol(graw))

gm=snp.df$map
names(gm)=snp.df$markerName
gmap.s=split(gm, snp.df$chr)
cross=buildCrossObject(graw,gmap.s,'f2','_')
gps=calc_genoprob(cross, error_prob=.001)

x11()
ind=5
gpi=do.call('cbind', lapply(gps, function(x) x[ind,,]))

par(mfrow=c(2,1))
plot(-1+real[,ind]*2, type='h', ylab='simulated genotype')
points(-1+2*graw[,ind], pch=20,cex=2, col='yellow')
#points((2*gpi[3,]),col='blue', ylab='p(0)')
#points((2*(-1*gpi[1,])),col='red', ylab='p(2)')
#points(2*gpi[2,]-1,col='grey', ylab='p(1)')
abline(v=cumsum(sapply(gmap.s, length)),lty=2)

plot(gpi[3,],col='blue', ylim=c(0,1), ylab='genotype probability')
points(gpi[2,],col='grey')
points(gpi[1,],col='red')
points(graw[,ind],cex=2, pch=20, col='yellow' )
abline(v=cumsum(sapply(gmap.s, length)),lty=2)


