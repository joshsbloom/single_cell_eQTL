 # Eskin trick to compute heritability faster ------------------------------------------------------------
    #helper functions
    doEigenA_forMM=function(pheno.scaled,A ,X=NULL ) {
            n=nrow(pheno.scaled)
            if(is.null(X) ) {  X = matrix(rep(1, n), n, 1); p=1 } else {p=ncol(X) }
            XtX = crossprod(X, X)
            XtXinv = solve(XtX)
            S = diag(n) - tcrossprod(X %*% XtXinv, X)
            SHbS = S %*% A %*% S
            SHbS.system = eigen(SHbS, symmetric = TRUE)
            theta = SHbS.system$values[1:(n - p)] 
            Q = SHbS.system$vectors[, 1:(n - p)]
            return(list(theta=theta, Q=Q))
            }
    # borrowed and modified from  from mixed.solve() in rrBLUP 08/02/15 -----------------------------------
    # super optimized for one VC and fixed effects ~1000X speedup by precomputing eigen decomp
    m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09), theta=NULL, Q=NULL, X=NULL ) 
    {
        n <- length(y)
        y <- matrix(y, n, 1)
        if(is.null(X) ) {  p <- 1    } else { p = ncol(X) }
        Z <- diag(n)
        m <- ncol(Z)
           
        omega <- crossprod(Q, y)
        omega.sq <- omega^2
        
        f.REML <- function(lambda, n.p, theta, omega.sq) {
            n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
        }
        soln <- optimize(f.REML, interval = bounds, n - p, theta,  omega.sq)
        lambda.opt <- soln$minimum
        
        df <- n - p
        Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
        Ve.opt <- lambda.opt * Vu.opt
        VCs=c(Vu.opt, Ve.opt)
        return(VCs)
    }
    calcA=function(p,A,do.print=T) {
        vcA.e=cbind(rep(NA, ncol(p)), rep(NA, ncol(p)))
        rownames(vcA.e)=colnames(p)
        eigA.e=doEigenA_forMM(p,A)
        # calculate mixed model, one term for additive variance  -------------------------------------------
        vcA.e=foreach(i=1:ncol(p), .combine='rbind') %dopar% {
                if(is.na(sd(p[,i]))) {
                   return(c(NA,NA))
             # next;
              }
            #vcA.e[i,]=m.S(p[,i], K=A,  theta=eigA.e$theta, Q=eigA.e$Q)
            return(m.S(p[,i], K=A,  theta=eigA.e$theta, Q=eigA.e$Q))
            #if(do.print){     print(i)}
        }
        return(vcA.e)
    }
#-------------------------------------------------------------------------------------------------------------
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%crossprod(Z,Vinv)%*%(y- X%*%B)    }

