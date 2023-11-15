library(MASS)
library(glmmTMB)
library(foreach)
library(doMC)
library(S4Vectors)
library(ggplot2)

registerDoMC(cores=64)

x=readRDS('/data/single_cell_eQTL/yeast/results/18_3051_May10/nbin_A.RDS')
m=x$estimate[x$term=='genoB' & x$effect=='fixed' & x$component=='cond']
d=x$estimate[x$term=='genoB' & x$effect=='fixed' & x$component=='disp']
p=x$dispLRT[x$term=='genoB' & x$effect=='fixed' & x$component=='disp']
# fig 4 
plot(log(exp(m[p<.001])),log(1/exp(d[p<.001])))


#plot(m[p<.001],d[p<.001])
#median allelic mean
quantile(attr(x, 'gm')$emmean, .5)
quantile(attr(x, 'gd')$emmean, .5)

#fold changes
hist(na.omit(x$estimate[x$term=='genoB' & x$effect=='fixed' & x$component=='cond' & x$dispLRTp<.01]))
hist(na.omit(x$estimate[x$term=='genoB' & x$effect=='fixed' & x$component=='disp' & x$dispLRTp<.01]))

# ~.05 counts per cell

#nseg.range=c(250,500,750,1000,2000)
nseg.range=c(5000)
#BCV.range=seq(.1,1,.2)
#count.range=2^seq(-5,2,1)
#roughly mean counts per transcript per cell and 95% quantile
mu.int=c(-3) #,-.7) #(c(.05,.5))
#theta.range=1/BCV.range^2 
fc.range=seq(-3,3,.5)  #  #c(-1.5, 1.5))

theta.int=c(-2)
theta.fc=seq(-3,3,.5)

sim.space=expand.grid(mu.int, fc.range, theta.int, theta.fc, nseg.range)
names(sim.space)=c('mu', 'fc', 'theta1', 'theta2', 'nsegs')
sim.space$p.diff=NA
sim.space$p.disp=NA

sim.pow=foreach(simIter=1:nrow(sim.space)) %dopar% { #nrow(sim.space)) %dopar% {
    print(paste(simIter, sim.space[simIter,]))
    ncells=sim.space[simIter, 'nsegs']
    mu=sim.space[simIter, 'mu']
    fc=sim.space[simIter, 'fc']
    theta1=sim.space[simIter, 'theta1']
    theta2=sim.space[simIter, 'theta2']
   # total.count.per.cell=sim.space[simIter, 'count']
    nrep=250
    p.diff=rep(NA,nrep)
    p.disp=p.diff #rep(NA,100)
    b.mean.int=p.diff #rep(NA,100)
    b.disp.int=p.diff #rep(NA,100)
    b.mean=p.diff #rep(NA,)
    b.disp=p.diff #rep(NA,
    for(i in 1:nrep){
        #print(i)
        ycounts1=rnegbin(ncells, exp(mu), theta=exp(theta1))
        ycounts2=rnegbin(ncells, exp(mu+fc),theta=exp(theta1+theta2))
        geno=factor(c(rep('A',ncells),rep('B',ncells)))
        cellID=factor(as.character(seq(1:ncells)))
        cellIDf=c(cellID,cellID)
        y=c(ycounts1,ycounts2)
        test=glmmTMB(y~geno,family=nbinom2, dispformula=~geno) 
        b.mean.int[i]=as.numeric(fixef(test)$cond)[1]
        b.mean[i]=as.numeric(fixef(test)$cond)[2]
        b.disp.int[i]=as.numeric(fixef(test)$disp)[1]
        b.disp[i]=as.numeric(fixef(test)$disp)[2]
        stc=summary(test)$coefficients
        p.disp[i]=stc$disp[2,4]
        p.diff[i]=stc$cond[2,4]
    }
    return(list(p.diff=p.diff, p.disp=p.disp, b.mean.int=b.mean.int, b.disp.int=b.disp.int, b.mean=b.mean, b.disp=b.disp))
}


pow.diff=sapply(sim.pow, function(x) {
       p.diff=x$p.diff
       sum(p.diff<.001,na.rm=T)/sum(!is.na(p.diff))
})
sim.space$p.diff=pow.diff
pow.disp=sapply(sim.pow, function(x) {
       p.disp=x$p.disp
       sum(p.disp<.001,na.rm=T)/sum(!is.na(p.disp))
})
sim.space$p.disp=pow.disp

sim.space$delta.mean=sapply(sim.pow, function(x) mean(x$b.mean))
sim.space$delta.disp=sapply(sim.pow, function(x) mean(x$b.disp))






ggplot(sim.space,
       aes(x=fc,
           y=delta.mean, 
           color=(delta.disp),
           group=(delta.disp)
           ))+
      labs(color='delta.disp')+ylab(c('delta.mean estimate'))+xlab('simulated mean delta')+
    geom_point()+geom_abline(intercept=0,slope=1)


ggplot(sim.space,
       aes(x=-theta2,
           y=-delta.disp, 
           color=factor(fc),
           group=fc
           ))+
      labs(color='mean fc')+ylab(c('delta.disp estimate'))+xlab('simulated disp delta')+
    geom_point()+geom_line()+geom_abline(intercept=0,slope=1)



ggplot(sim.space,
       aes(x=fc,
           y=-theta2, 
           color=factor(fc),
           group=fc
           ))+
      labs(color='mean fc')+ylab(c('simulated disp delta'))+xlab('simulated mean delta')+
    geom_point()


ggplot(sim.space,
       aes(x=delta.mean,
           y=-delta.disp, 
           color=factor(fc),
           group=fc
           ))+
      labs(color='mean fc')+ylab(c('delta.disp estimate'))+xlab('delta.mean estimate')+
    geom_point()

geom_abline(intercept=0,slope=1)

jsim.means=unlist(lapply(sim.pow, function(x) x$b.mean))
jsim.disps=unlist(lapply(sim.pow, function(x) x$b.disp))
jsim.dispp=unlist(lapply(sim.pow, function(x) x$p.disp))


