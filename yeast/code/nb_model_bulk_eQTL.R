library(glmmTMB)
library(lme4)
library(MASS)
library(edgeR)
library(regress)
library(foreach)
library(Rfast)
library(doMC)
registerDoMC(cores=64)

source('/data/single_cell_eQTL/yeast/code/helper_h2_functions.R')

# important objects from bulk eQTL experiment
# saveRDS(count.matrix, file='/data/single_cell_eQTL/yeast/Bulk/data/countMatrix.RDS')
# saveRDS(peaks.per.gene, file='/data/single_cell_eQTL/yeast/Bulk/data/ppg.RDS')
# saveRDS(gbatch.fact, file='/data/single_cell_eQTL/yeast/Bulk/data/batch_factor.RDS')
# saveRDS(OD.cov, file='/data/single_cell_eQTL/yeast/Bulk/data/OD_cov.RDS')
# saveRDS(gdata, file='/data/single_cell_eQTL/yeast/Bulk/data/gdata.RDS')
gdata=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/gdata.RDS')
count.matrix=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/countMatrix.RDS')
peaks.per.gene=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/ppg.RDS')
gbatch.fact=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/batch_factor.RDS')
OD.cov=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/OD_cov.RDS')
#library(rio)
#export(gdata, '/data/single_cell_eQTL/yeast/Bulk/data/genotypes.csv', row.names=T)
#export(count.matrix, '/data/single_cell_eQTL/yeast/Bulk/data/countMatrix.csv',row.names=T)
#export(gbatch.fact, '/data/single_cell_eQTL/yeast/Bulk/data/batch_factor.csv', row.names=T)
#technical_covariates=data.frame(model.matrix(gdata[,1]~gbatch.fact+OD.cov-1))
#export(technical_covariates, '/data/single_cell_eQTL/yeast/Bulk/data/technical_covariates.csv')

Amat=tcrossprod(scale(gdata))/(ncol(gdata)-1)


CPM=cpm(count.matrix, normalized.lib.sizes=T)
aveLogCPM=apply(CPM, 1, function(x) mean(log2(x+1)))


tpm.matrix=log2(apply(count.matrix,2, function(x) countToTpm(x, gene.annot.df$length))+0.5)


# residuals from count-based model for heritability estimate
count.resids=matrix(NA, nrow(count.matrix), ncol(count.matrix))
colnames(count.resids)=colnames(count.matrix)
rownames(count.resids)=rownames(count.matrix)

#augment peak data with effect size estimates from count based regression model 
peaks.per.geneGLM=list()
for(gene in names(peaks.per.gene)){
    print(gene)
    #gene=names(peaks.per.gene)[11]#'YJR009C'
    y=count.matrix[gene,]
    strain.counts=colSums(count.matrix)
    XX=gdata[,peaks.per.gene[[gene]]$pmarker]
    
    # deal with glm.nb convergence fail by fitting poisson regression instead
    #m0=lm(log2(y)~offset(log(strain.counts))+gbatch.fact+scale(OD.cov)+XX)
    #m1=tryCatch({
    #     glm.nb(y~offset(log(strain.counts))+gbatch.fact+scale(OD.cov)+XX) },
    #     error=function(e) {
    #        glm(y~offset(log(strain.counts))+gbatch.fact+scale(OD.cov)+XX,family='poisson')
    #    })
    m2=residuals(glm.nb(y~offset(log(strain.counts))+gbatch.fact+scale(OD.cov)))
    count.resids[gene,]=as.vector(m2)
   # tm=tidy(m1)[-c(1:14),]
   # tm$term=gsub('^XX', '', tm$term)
    # peaks.per.geneGLM[[gene]]=cbind(peaks.per.gene[[gene]],tm)
}

# should be good enough for now
h2_qtl=calcA(t(count.resids), Amat)
rownames(h2_qtl)=rownames(count.resids)
colnames(h2_qtl)=c('A', 'E')

edg






A2=data.matrix(nearPD(Amat,corr=T)$mat)
A2.inv=solve(A2)
strain.counts=colSums(count.matrix)
off=log(strain.counts)
Zid=as.factor(colnames(A2))
data2=data.frame(off=off, gbatch.fact=gbatch.fact, OD.cov=OD.cov, Zid=Zid)


x=readRDS('/data/single_cell_eQTL/yeast/results/combined/Ap/scICC_3VC.RDS')


 ICC.R=sR/(sA+sC+sR+log(1+(1/lmbda)+(1/sgma)))
        ICC.C=sC/(sA+sC+sR+log(1+(1/lmbda)+(1/sgma)))


bGLMM_bulk=list()
    for(gene in rownames(x)[rownames(x) %in% names(peaks.per.gene)] ) { # names(peaks.per.gene)){
        print(gene)
        y=round(count.matrix[gene,]  )
        names(y)=colnames(A2)
        bGLMM_bulk[[gene]]=spaMM::fitme(y~offset(off)+gbatch.fact+OD.cov+corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2)
        
        print(x[gene,])
        print(VarCorr( bGLMM_bulk[[gene]] ))
    }

saveRDS(bGLMM_bulk, file='/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMM_bulk.RDS')

lab=rep(0, length(bGLMM_bulk))
for(i in 1:length(bGLMM_bulk)) {
    print(i)
    f3=bGLMM_bulk[[i]]
    #lmbda=mean(exp(predict(f3)))
    lab[i]=mean(epf[is.finite(epf)])
    }



vab=sapply(bGLMM_bulk, function(x) VarCorr(x)[1,3])

saveRDS(vab, file='/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMM_bulk_A.RDS')

vab=readRDS('/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMM_bulk_A.RDS')

dab=rep(NA,length(bGLMM_bulk))
for( i in 1:length(bGLMM_bulk) ){
    test=bGLMM_bulk[[i]]
    a=capture.output(summary(test))[4]
    dab[i]=as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", a))
    print(dab[i])
}
names(dab)=names(vab)
saveRDS(dab, file='/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMM_bulk_theta.RDS')

dab=readRDS('/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMM_bulk_theta.RDS')


veb=log(1+(1/lab)+(1/dab))


nb_bulk_res=data.frame(gene=names(vab), Var.A=vab, Var.E=veb, lmbda=lab, sgma=dab, h2=vab/(vab+veb))

sub=vcA.OD.unscaled[match(nb_bulk_res$gene, rownames(vcA.OD.unscaled)),]



dab =sapply(bGLMM_bulk, function(test) { as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", capture.output(summary(test)))[3]) })


for(i in 1:length(x)){

}


xab=rep(NA,length(x))
for( i in 1:length(x) ){
    test=x[[i]]
    a=capture.output(summary(test))[5]
    xab[i]=as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", a))
    print(xab[i])
}
names(xab)=names(x)
saveRDS(xab, file='/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMMsA_theta.RDS')

xab=rep(NA,length(x))
for( i in 1:length(x) ){
    test=x[[i]]
    a=capture.output(summary(test))[5]
    xab[i]=as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", a))
    print(xab[i])
}
names(xab)=names(x)
saveRDS(xab, file='/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMMsA_theta.RDS')

xab=readRDS('/data/single_cell_eQTL/yeast/results/combined/Ap/bGLMMsA_theta.RDS')







#f= #Poisson(link="log"))
        print(colnames(Yr)[i])
        bGLMMsA[[colnames(Yr)[i]]]=
         
        print(bGLMMsA[[colnames(Yr)[i]]])




library(MASS)
library(mvnfast)
library(spaMM)

y=rnegbin(1012,exp(m[1,])*10,theta=100)
sigma(glmmTMB(y~1,family=nbinom2))
rcov=.4*A2
m=rmvn(1,mu=rep(0,1012), rcov)

y=rnegbin(1012,exp(m[1,]),theta=100)
system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2())})
system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=F))})
system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=F, NbThreads=36))})
system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=T, NbThreads=36))})

system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(corrMatrix=A2), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=T, NbThreads=72))})
system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=T, NbThreads=72))})

system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(corrMatrix=A2), method='PQL', data=data2, family=negbin2(), control.HLfit=list(sparse_precision=F, NbThreads=36))})





mur=c(.01, 0.1, 1, 10, 100, 100)
thetar=c(.5, 1, 5, 10, 50, 100)

eg=expand.grid(mur,thetar)

biasL=list()
for(i in 1:nrow(eg)){
    print(i)
biasL[[as.character(i)]]=replicate(100, {
                                       
            y=rnegbin(1012,rep(eg[i,1], 1012) ,theta=eg[i,2])
            test=spaMM::fitme(y~1, method='PQL', family=negbin2(),data=data.frame(y=y)) 
            a=capture.output(summary(test))[4]
            return(as.numeric(gsub(".*shape=(.*)\\)\\( link.*", "\\1", a)))
            }) 
print(biasL[[as.character(i)]])
    
}


test=











system.time({test=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='ML', data=data2, family=negbin2)})
system.time({test2=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='REML', data=data2, family=negbin2)})
system.time({test3=spaMM::fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL/L', data=data2, family=negbin2)})

system.time({test3=fitme(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='ML', data=data2, family=negbin())})

glm.nb(y~gdata[,seq(1,10000,100)])$theta
sigma(glmmTMB(y~gdata[,seq(1,10000,100)],family=nbinom2))

system.time({test3=corrHLfit(y~corrMatrix(1|Zid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin)})

ly=log(y+1)
r=regress(ly~1,~A2,verbose=F)

  lmbda=mean(exp(predict(test)))
        sgma=sigma(test)
        sA=VarCorr(test)[1,3]
       # sR=VarCorr(test)[2,3]
       # sC=VarCorr(test)[3,3]
r$sigma[1]/(sum(r$sigma))

sA/(sA+log(1+(1/lmbda)+(1/sgma)))
sA/(sA+trigamma(1/((1/lmbda)+(1/sgma))))
 #trigamma(1/((1/lmbda)+(1/sgma)))) #log(1+(1/lmbda)+(1/sgma)))
































#saveRDS(peaks.per.geneGLM, file='/data/single_cell_eQTL/yeast/Bulk/data/peaks.per.geneGLM.RDS')
#saveRDS(h2_qtl, file='/data/single_cell_eQTL/yeast/Bulk/data/h2_qtl.RDS')
#saveRDS(count.resids, file='/data/single_cell_eQTL/yeast/Bulk/data/countResids.RDS')

# refit and estimate betas given a count-based model
# extract Betas
eqtl_betas=list()
eqtl_h2=list()
for(gene in names(peaks.per.gene)){
    print(gene)
    
    sg=scale(gdata[,peaks.per.gene[[gene]]$pmarker])
    gm=lm(scale(count.resids[gene,])~sg)
    h2=summary(gm)$r.squared
    cmm=coef(gm)[-1]
    names(cmm)=peaks.per.gene[[gene]]$pmarker #gsub('^sg','', names(cmm))

    eqtl_betas[[gene]]=cmm
    eqtl_h2[[gene]]=h2
}

saveRDS(eqtl_betas, file='/data/single_cell_eQTL/yeast/Bulk/data/eqtl_betas.RDS')
saveRDS(eqtl_h2, file='/data/single_cell_eQTL/yeast/Bulk/data/eqtl_h2.RDS')

