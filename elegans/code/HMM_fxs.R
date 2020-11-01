 # for F2
mTmat=function(r){    matrix(c((1-r)^2,2*r*(1-r),r^2,r*(1-r),(1-r)^2+r^2, r*(1-r), r^2, 2*r*(1-r),(1-r)^2),3,3)  }
# for haploids
mTmatH=function(r) {   matrix(c((1-r),r,r,(1-r)),2,2) }

#convert haldane genetic map dist back to recombination fraction
mf.h = function(d) 0.5*(1-exp(-d/50))
map2rf = function(map, tol=1e-12) {
    rf = mf.h(diff(map))
    rf[rf < tol] = tol # don't let values be too small
    rf
}

# forward algorithm
calcForward=function(eMats, tMats, startProbs=NULL,lik=F) {
    if(is.null(startProbs)){
        if(nrow(eMats)==3){ startProbs=c(.25,.5,.25) }
        if(nrow(eMats)==2){ startProbs=c(.5,.5)      }
    }
    emissionProbs=eMats
    nObservations=ncol(eMats)
    nStates=nrow(eMats)
    f=array(NA, c(nStates, nObservations))
    States=seq(1:nrow(eMats))
    dimnames(f)=list(states=rownames(eMats), index=colnames(eMats))         
    logsumTot=0
   
    for (state in States) {     f[state, 1] = sum(log(startProbs[state] * emissionProbs[state, 1]))     }
    for (kk in 2:nObservations) {
            for (state in States) {
                logsum = -Inf
                for (previousState in States) {
                    temp = f[previousState, kk - 1] + log(tMats[[kk-1]][previousState,state]) 
                    if (temp > -Inf) {
                        if(!is.finite(exp(logsum-temp))) {
                            logsum = logsum
                        } else {
                            logsum = temp + log(1 + exp(logsum - temp))
                        }
                    }
                }
                f[state, kk] = sum(log(emissionProbs[state, kk])) + logsum
                if(lik) { if(is.finite(f[state, kk])) { logsumTot=logsumTot + f[state, kk] } }
            }
    }
    if(lik){return(logsumTot)}
    #print(logsumTot)
    return(f)
}
# backward algorithm
calcBackward=function(eMats, tMats) {
    emissionProbs=eMats
    nObservations=ncol(eMats)
    nStates=nrow(eMats)
    States=seq(1:nrow(eMats))

    b=array(NA, c(nStates, nObservations))
    dimnames(b)=list(states=rownames(eMats), index=colnames(eMats)) 

    for (state in States) {   b[state, nObservations] = log(1)    }
    for (k in (nObservations - 1):1) {
            for (state in States) {
                logsum = -Inf
                for (nextState in States) {
                    temp = b[nextState, k + 1] + sum(log(tMats[[k]][state, nextState] * emissionProbs[nextState, k + 1]))
                    if (temp > -Inf) {
                        if(!is.finite(exp(logsum-temp))) {
                          logsum = logsum
                         } else {  
                             logsum = temp + log(1 + exp(logsum - temp)) 
                        }
                    }
                }
                b[state, k] = logsum
            }
    }
    return(b)
}

#viterbi path
viterbi=function (eMats, tMats, startProbs=NULL) {
     if(is.null(startProbs)){
        if(nrow(eMats)==3){ startProbs=c(.25,.5,.25) }
        if(nrow(eMats)==2){ startProbs=c(.5,.5)      }
    }
    emissionProbs=eMats
    nObservations=ncol(eMats)
    nStates=nrow(eMats)
    States=seq(1:nrow(eMats))
    v=array(NA, c(nStates, nObservations))
    dimnames(v)=list(states=rownames(eMats), index=colnames(eMats)) 

    for (state in States) {
        v[state, 1] = sum(log(startProbs[state] * emissionProbs[state, 1]))
    }
    for (kk in 2:nObservations) {
        for (state in States) {
            maxi = NULL
            for (previousState in States) {
                temp = v[previousState, kk - 1] + log(tMats[[kk-1]][previousState,state]) 
                maxi = max(maxi, temp)
            }
            v[state, kk] = sum(log(emissionProbs[state, kk])) +  maxi
        }
    }
    viterbiPath = rep(NA, nObservations)
    for (state in States) {
        if (max(v[, nObservations]) == v[state, nObservations]) {
            viterbiPath[nObservations] = state
            break
        }
    }
    for (kk in (nObservations - 1):1) {
        for (state in States) {
            if (max(v[, kk] + log(tMats[[kk]][, viterbiPath[kk + 1]])) == v[state, kk] + log(tMats[[kk]][state,  viterbiPath[kk + 1]])) {
                viterbiPath[kk] = state
                break
            }
        }
    }
    return(viterbiPath)
}


#calculate posterior probabilities
calcPosterior=function(f,b) {
    nObservations=ncol(f)
    nStates=nrow(f)
    probObservations = f[1, nObservations]
    for (i in 2:nStates) {
        j = f[i, nObservations]
        if (j > -Inf) {
            probObservations = j + log(1 + exp(probObservations - j))
        }
    }
    posteriorProb = exp((f + b) - probObservations)
    return(posteriorProb)
}


