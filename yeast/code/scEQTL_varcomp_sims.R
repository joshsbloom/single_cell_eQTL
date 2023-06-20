 #custom sims
A=tcrossprod(scale(prevG))/ncol(prevG)
cov=.03*A+.000001*diag(nrow(A)) # #+.2*diag(nrow(Z))
test.theta=8
mvnormsim=mvnfast::rmvn(10000, mu=rep(0, nrow(cov)), sigma=cov)
mu=t(apply(mvnormsim ,1 , exp))
Y=apply(mu, 1, function(x) MASS::rnegbin(nrow(cov), mu=x, theta=test.theta))
Ysum=rowSums(Y)
library(regress)

lYsum=scale(log2(Ysum+1))
regress(lYsum~1, ~A, verbose=T)


thetas=sapply(bGLMMs, sigma)
thetas[thetas>1000]=1000
Z=model.matrix(Yr[,1]~Zid-1)
C=model.matrix(Yr[,1]~Cid-1)
ZtZ=Z%*%t(Z)
CtZ=C%*%t(C)

sus=sort(unique(best_match_seg))
A=tcrossprod(scale(t(sprevG[,sus])))/nrow(sprevG)
A2=data.matrix(nearPD(A,corr=T)$mat)
A2.inv=solve(A2)
ZAtZ=Z%*%A%*%t(Z)

test.theta=8
sAv=rep(NA,100)
sCv=rep(NA,100)
for(i in 1:100) {
    print(i)
    cov=.2*ZtZ+.1*CtZ+.4*ZAtZ+.00001*diag(nrow(Z)) # #+.2*diag(nrow(Z))
    mvnormsim=mvnfast::rmvn(1, mu=rep(0, nrow(cov)), sigma=cov)
    y=MASS::rnegbin(nrow(cov), mu=exp(mvnormsim[1,]), theta=test.theta)
    system.time({test=glmmTMB(y~(1|Cid)+(1|Zid), family=nbinom2, control=glmmTMBControl(parallel=36))})
    summary(test)
    sA=VarCorr(test)$cond$Zid[1]
    print(sA)
    sAv[i]=sA
    sC=VarCorr(test)$cond$Cid[1]
    print(sC)
    sCv[i]=sC
}
lmbda=mean(exp((predict(test))))
sgma=sigma(test)
sC/(sA+sC+log(1+(1/lmbda)+(1/sgma)))
sA/(sA+sC+log(1+(1/lmbda)+(1/sgma)))

test2=glmer.nb(y~(1|Cid)+(1|Zid))
data2=data
data2$y=y

#f=fitme(y~(1|Zid)+ corrMatrix(1|Zid), corrMatrix=A2, method='REML', family=Poisson(link='log'),data=data2)
#f=fitme(y~1+(1|Zid)+corrMatrix(1|Zid), covStruct=list(corrMatrix=A2), method='REML', data=data2, family=Poisson(link="log"))
system.time({f=fitme(y~1+corrMatrix(1|Zid)+(1|Zid)+(1|Cid), covStruct=list(precision=A2.inv), method='PQL', data=data2, family=negbin2)}) #Poisson(link="log"))})

as.data.frame(VarCorr(test2))
