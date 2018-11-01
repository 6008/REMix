MCEM <- function(Int,Mu,Sigma,V,Pi,verb=TRUE,niter,probs,q,knots=knots){
  ## i: read
  ## j: cycle
  ## k: channel (A,C,G,T)
  library(mvtnorm)
  library(abind)
  library(Matrix)
  nread <- dim(Int)[1]
  nchannel <- dim(Int)[2]
  ncycle <- dim(Int)[3]
  ## creating the design matrices for each cycle j
  dimX <- length(knots)*4+2*4
  pred <- 0:(ncycle-1)
  kfn <- list()
  for(kn in 1:length(knots)){
    kfn[[kn]] <- ifelse(pred>knots[kn],pred-knots[kn],0)
  }
  Xj <- list()
  for(j in 1:ncycle){
    tmp.X1 <- cbind(diag(1,nchannel),diag(j-1,nchannel))
    tmp.X2 <- NULL
    for(kn in 1:length(knots)){
      tmp.X2 <- cbind(tmp.X2,diag(kfn[[kn]][j],nchannel))
    }
    Xj[[j]] <- cbind(tmp.X1,tmp.X2)
  }

  bigX <- Xj[[1]]
  for (j in 2:ncycle){
    bigX <- rbind(bigX, Xj[[j]])
  }

  PiAvg <- rep(0,nchannel)
  MuAvg <- matrix(0,ncol=4,nrow=dimX)
  SigmaAvg <- list()
  SigmaAvg[[1]] <- SigmaAvg[[2]] <- SigmaAvg[[3]] <- SigmaAvg[[4]] <- matrix(0,ncol=dimX,nrow=dimX)
  VAvg <- list()
  VAvg[[1]] <- VAvg[[2]] <- VAvg[[3]] <- VAvg[[4]] <- matrix(0,ncol=4,nrow=4)
  l <- 0

  P <- probs
  
  iter <- 0

      ## Generating the deltas
    ptm <- proc.time()

    Delta <- array(0,dim=c(nread,nchannel,ncycle))
    for(i in 1:nread){
      for(j in 1:ncycle){
        Delta[i,,j] <- t(rmultinom(1,1,P[i,,j]))
      }
    }

    print("Generating the deltas: Done")
    print(proc.time() - ptm)

    ## Getting the indices for the cycles that belong to each read and component
    ptm <- proc.time()

    IndByGroup <- list()
    for(k in 1:nchannel){
      IndByGroup[[k]] <- list()
      for(i in 1:nread){
        IndByGroup[[k]][[i]] <- list()
        IndByGroup[[k]][[i]] <- which(Delta[i,k,]==1)
      }
    }

    print("Getting the indices for the cycles that belong to each read and component: Done")
    print(proc.time() - ptm)

    ## obtaining the covariance matrices for the conditional distributions
    ptm <- proc.time()

    BigV <- list()
    BigV[[1]] <- V[[1]]
    BigV[[2]] <- V[[2]]
    BigV[[3]] <- V[[3]]
    BigV[[4]] <- V[[4]]
    for(j in 2:ncycle){
      BigV[[1]] <- bdiag(BigV[[1]],V[[1]])
      BigV[[2]] <- bdiag(BigV[[2]],V[[2]])
      BigV[[3]] <- bdiag(BigV[[3]],V[[3]])
      BigV[[4]] <- bdiag(BigV[[4]],V[[4]])
    }
    BigSigma12 <- BigSigma22 <- BigSigma22Inverse <- BigSigma22Rows <- list()
    for(k in 1:nchannel){
      BigSigma12[[k]] <- list()
      BigSigma22Rows[[k]] <- list()
      for(j in 1:ncycle){
        BigSigma12[[k]] <- cbind(BigSigma12[[k]],Sigma[[k]]%*%t(Xj[[j]]))
      }
      for(j in 1:ncycle){
        BigSigma22Rows[[k]][[j]] <- list()
        BigSigma22Rows[[k]][[j]] <- Xj[[j]]%*%matrix(as.numeric(BigSigma12[[k]]),ncol=(nchannel*ncycle), byrow=FALSE)
      }
      BigSigma22[[k]] <- abind(BigSigma22Rows[[k]],along=1) + as.matrix(BigV[[k]])
      BigSigma12[[k]] <- matrix(as.numeric(BigSigma12[[k]]),ncol=(nchannel*ncycle),byrow=FALSE)
    }

    print("Obtaining the covariance matrices for the conditional distributions: Done")
    print(proc.time() - ptm)

    ## Calculating the posterior means
    ptm <- proc.time()
    E.gamma <- list()
    Cov.gamma <- list()
    Sigma12 <- list()
    S12S22 <- list()
    for(k in 1:nchannel){
      E.gamma[[k]] <- matrix(0,ncol=nread,nrow=dimX)
      Cov.gamma[[k]] <- list()
      for(i in 1:nread){
        if(length(IndByGroup[[k]][[i]])>0){
          CycleIndexi <- c(c(IndByGroup[[k]][[i]])*4,c(IndByGroup[[k]][[i]])*4-3,c(IndByGroup[[k]][[i]])*4-2,c(IndByGroup[[k]][[i]])*4-1)
          Sigma12[[k]] <- BigSigma12[[k]][,CycleIndexi]
          S12S22[[k]] <- Sigma12[[k]] %*% solve(BigSigma22[[k]][CycleIndexi,CycleIndexi])
          Cov.gamma[[k]][[i]] <- matrix(0, ncol=dimX, nrow=dimX)
          E.gamma[[k]][,i] <- Mu[,k] + S12S22[[k]] %*% (c(Int[i,,]) - bigX %*% Mu[,k])[CycleIndexi,1]
          Cov.gamma[[k]][[i]] <- Sigma[[k]] - S12S22[[k]] %*% t(Sigma12[[k]])
        }else{
          E.gamma[[k]][,i] <- Mu[,k]
          Cov.gamma[[k]][[i]] <- Sigma[[k]]
        }
      }
    }
    print("Calculating the posterior means: Done")
    print(proc.time() - ptm)

    ## Step 2: Update mixture proportions
    ptm <- proc.time()

    Pi.new <- apply(Delta,2,sum) / (nread*ncycle)

    print("Step 2 (Pi): Done")
    print(proc.time() - ptm)

    ## Step 3: Update mean of random effects
    ptm <- proc.time()

    Mu.new <- matrix(0,ncol=nchannel,nrow=dimX)
    for(k in 1:nchannel){
      Mu.new[,k] <- apply(E.gamma[[k]],1,sum)
    }
    Mu.new <- Mu.new / nread

    print("Step 3 (Mu): Done")
    print(proc.time() - ptm)

    ## Step 4: Update covariance of random effects
    ptm <- proc.time()

    Sigma.new <- list()
    for(k in 1:nchannel){
      Sigma.new[[k]] <- (Reduce("+",Cov.gamma[[k]]) + (E.gamma[[k]]-Mu.new[,k]) %*% t(E.gamma[[k]]-Mu.new[,k])) / nread
    }

    print("Step 4 (Sigma): Done")
    print(proc.time() - ptm)

    ## Step 5: Update residual variance
    ptm <- proc.time()

    V.new <- list()
    V.new[[1]] <- V.new[[2]] <- V.new[[3]] <- V.new[[4]] <- matrix(0,nrow=4,ncol=4)
    ## j 1:ncycles
    for(k in 1:nchannel){
      for(j in 1:nchannel){
        E.temp <- t(Int[which(Delta[,k,j]==1),,j]) - Xj[[j]] %*% E.gamma[[k]][,which(Delta[,k,j]==1)]
        Cov.temp <- as.matrix(Reduce("+",Cov.gamma[[k]][which(Delta[,k,j]==1)]))
        V.new[[k]] <- V.new[[k]] + (Xj[[j]] %*% Cov.temp %*% t(Xj[[j]]) + E.temp %*% t(E.temp))
      }
      V.new[[k]] <- V.new[[k]] / sum(Delta[,k,])
    }

    print("Step 5 (ResVar): Done")
    print(proc.time() - ptm)

    iter <- iter + 1
    Pi <- Pi.new
    Mu <- Mu.new
    Sigma <- Sigma.new
    V <- V.new

    if(verb){
      cat(paste("Iteration= ",iter,"\n",sep=""))
    }


  ############################################################################################################
  ############################################################################################################
  while (iter < niter){
    ## Step 1.a (E-Step): Generate Deltas
    ## Calcuting mixture probabilities
    ptm <- proc.time()

    ##pi is really pi
    P.denom <- array(0,dim=c(nread,nchannel,ncycle))
    for(k in 1:nchannel){
      for(j in 1:ncycle){
        Mean <- Xj[[j]]%*%Mu[,k]
        denom <- (2*pi)^(-nchannel/2) * det( V[[k]] + Xj[[j]] %*% Sigma[[k]] %*% t(Xj[[j]]) )^(-1/2)
        SigmaInv <- solve( V[[k]] + Xj[[j]] %*% Sigma[[k]] %*% t(Xj[[j]]) )
        for(i in 1:nread){
          P.denom[i,k,j] <- denom * exp(-0.5*t(Int[i,,j]-Mean) %*% SigmaInv %*% (Int[i,,j]-Mean) )
        }
      }
    }
    P1.tmp <- array(0,dim=c(nread,nchannel,ncycle))
    for(i in 1:nread){  
      for(j in 1:ncycle){
        P1.tmp[i,,j] <- P.denom[i,,j]*Pi
      }
    }
    P1.tmp.sum <- apply(P1.tmp,c(1,3),sum)
    qq <- quantile(P1.tmp.sum,q) ## trimming data
    DensityInd <- which(P1.tmp.sum<qq,arr.ind=TRUE)
    CyclesForRead <- list()
    for(i in 1:nread){
      if(length(DensityInd[DensityInd[,1]==i,2])>0){
        CyclesForRead[[i]] <- (1:ncycle)[-DensityInd[DensityInd[,1]==i,2]]
      }else{
        CyclesForRead[[i]] <- 1:ncycle
      }
    }
    
    ## NOUSE?
    ReadsForCycle <- list()
    for(j in 1:ncycle){
      if(length(DensityInd[DensityInd[,2]==j,1])>0){
        ReadsForCycle[[j]] <- (1:nread)[-DensityInd[DensityInd[,2]==j,1]]
      }else{
        ReadsForCycle[[j]] <- 1:nread
      }
    }
    ## NOUSE?
    
    
    P <- array(0,dim=c(nread,nchannel,ncycle))
    for(i in 1:nread){
      for(j in CyclesForRead[[i]]){
        P[i,,j] <- P1.tmp[i,,j] / sum(P1.tmp[i,,j])
      }
    }

    print("Calcuting mixture probabilities: Done")
    print(proc.time() - ptm)

    ## Generating the deltas
    ptm <- proc.time()

    Delta <- array(0,dim=c(nread,nchannel,ncycle))
    for(i in 1:nread){
      for(j in CyclesForRead[[i]]){
        Delta[i,,j] <- t(rmultinom(1,1,P[i,,j]))
      }
    }

    print("Generating the deltas: Done")
    print(proc.time() - ptm)

    ## Getting the indices for the cycles that belong to each read and component
    ptm <- proc.time()

    IndByGroup <- list()
    for(k in 1:nchannel){
      IndByGroup[[k]] <- list()
      for(i in 1:nread){
        IndByGroup[[k]][[i]] <- list()
        IndByGroup[[k]][[i]] <- which(Delta[i,k,]==1)
      }
    }

    print("Getting the indices for the cycles that belong to each read and component: Done")
    print(proc.time() - ptm)

    ## obtaining the covariance matrices for the conditional distributions
    ptm <- proc.time()

    BigV <- list()
    BigV[[1]] <- V[[1]]
    BigV[[2]] <- V[[2]]
    BigV[[3]] <- V[[3]]
    BigV[[4]] <- V[[4]]
    for(j in 2:ncycle){
      BigV[[1]] <- bdiag(BigV[[1]],V[[1]])
      BigV[[2]] <- bdiag(BigV[[2]],V[[2]])
      BigV[[3]] <- bdiag(BigV[[3]],V[[3]])
      BigV[[4]] <- bdiag(BigV[[4]],V[[4]])
    }
    BigSigma12 <- BigSigma22 <- BigSigma22Inverse <- BigSigma22Rows <- list()
    for(k in 1:nchannel){
      BigSigma12[[k]] <- list()
      BigSigma22Rows[[k]] <- list()
      for(j in 1:ncycle){
        BigSigma12[[k]] <- cbind(BigSigma12[[k]],Sigma[[k]]%*%t(Xj[[j]]))
      }
      for(j in 1:ncycle){
        BigSigma22Rows[[k]][[j]] <- list()
        BigSigma22Rows[[k]][[j]] <- Xj[[j]]%*%matrix(as.numeric(BigSigma12[[k]]),ncol=(nchannel*ncycle), byrow=FALSE)
      }
      BigSigma22[[k]] <- abind(BigSigma22Rows[[k]],along=1) + as.matrix(BigV[[k]])
      BigSigma12[[k]] <- matrix(as.numeric(BigSigma12[[k]]),ncol=(nchannel*ncycle),byrow=FALSE)
    }
 
    print("Obtaining the covariance matrices for the conditional distributions: Done")
    print(proc.time() - ptm)

    ## Calculating the posterior means
    ptm <- proc.time()

    E.gamma <- list()
    Cov.gamma <- list()
    Sigma12 <- list()
    S12S22 <- list()
    for(k in 1:nchannel){
      E.gamma[[k]] <- matrix(0,ncol=nread,nrow=dimX)
      Cov.gamma[[k]] <- list()
      for(i in 1:nread){
        if(length(IndByGroup[[k]][[i]])>0){
          CycleIndexi <- c(c(IndByGroup[[k]][[i]])*4,c(IndByGroup[[k]][[i]])*4-3,c(IndByGroup[[k]][[i]])*4-2,c(IndByGroup[[k]][[i]])*4-1)
          Sigma12[[k]] <- BigSigma12[[k]][,CycleIndexi]
          S12S22[[k]] <- Sigma12[[k]] %*% solve(BigSigma22[[k]][CycleIndexi,CycleIndexi])
          Cov.gamma[[k]][[i]] <- matrix(0, ncol=dimX, nrow=dimX)
          E.gamma[[k]][,i] <- Mu[,k] + S12S22[[k]] %*% (c(Int[i,,]) - bigX %*% Mu[,k])[CycleIndexi,1]
          Cov.gamma[[k]][[i]] <- Sigma[[k]] - S12S22[[k]] %*% t(Sigma12[[k]])
        }else{
          E.gamma[[k]][,i] <- Mu[,k]
          Cov.gamma[[k]][[i]] <- Sigma[[k]]
        }
      }
    }
    print("Calculating the posterior means: Done")
    print(proc.time() - ptm)

    ## Step 2: Update mixture proportions
    ptm <- proc.time()

    Pi.new <- apply(Delta,2,sum) / sum(Delta)

    print("Step 2 (Pi): Done")
    print(proc.time() - ptm)

    ## Step 3: Update mean of random effects
    ptm <- proc.time()

    Mu.new <- matrix(0,ncol=nchannel,nrow=dimX)
    for(k in 1:nchannel){
      Mu.new[,k] <- apply(E.gamma[[k]],1,sum)
    }
    Mu.new <- Mu.new / nread

    print("Step 3 (Mu): Done")
    print(proc.time() - ptm)

    ## Step 4: Update covariance of random effects
    ptm <- proc.time()

    Sigma.new <- list()
    for(k in 1:nchannel){
      Sigma.new[[k]] <- (Reduce("+",Cov.gamma[[k]]) + (E.gamma[[k]]-Mu.new[,k]) %*% t(E.gamma[[k]]-Mu.new[,k])) / nread
    }

    print("Step 4 (Sigma): Done")
    print(proc.time() - ptm)

    ## Step 5: Update residual variance
    ptm <- proc.time()

    V.new <- list()
    V.new[[1]] <- V.new[[2]] <- V.new[[3]] <- V.new[[4]] <- matrix(0,nrow=4,ncol=4)
    for(k in 1:nchannel){
      for(j in 1:nchannel){
        E.temp <- t(Int[which(Delta[,k,j]==1),,j]) - Xj[[j]] %*% E.gamma[[k]][,which(Delta[,k,j]==1)]
        Cov.temp <- as.matrix(Reduce("+",Cov.gamma[[k]][which(Delta[,k,j]==1)]))
        V.new[[k]] <- V.new[[k]] + (Xj[[j]] %*% Cov.temp %*% t(Xj[[j]]) + E.temp %*% t(E.temp))
      }
      V.new[[k]] <- V.new[[k]] / sum(Delta[,k,])
    }

    print("Step 5 (V): Done")
    print(proc.time() - ptm)

    iter <- iter + 1
    Pi <- Pi.new
    Mu <- Mu.new
    Sigma <- Sigma.new
    V <- V.new
    
    if(verb){
      cat(paste("Iteration= ",iter,"\n",sep=""))
    }
    if (iter == niter){
      cat("Reached the maximum number of iterations", "\n")
    }
    if(iter>=(niter/2+1)){
      l <- l + 1
      PiAvg <- PiAvg + Pi
      MuAvg <- MuAvg + Mu
      for(k in 1:nchannel){
        SigmaAvg[[k]] <- SigmaAvg[[k]] + Sigma[[k]]
        VAvg[[k]] <- VAvg[[k]] + V[[k]]
      }
    }
  } ## ends while loop
  
  ## get final Pijk's for basecalling assignment
  Pi <- PiAvg / l
  Mu <- MuAvg / l
  Sigma <- list()
  V <- list()
  for(k in 1:nchannel){
    Sigma[[k]] <- SigmaAvg[[k]] / l
    V[[k]] <- VAvg[[k]] / l
  }
  Z.denom <- array(0,dim=c(nread,nchannel,ncycle))
  for(k in 1:nchannel){
    for(j in 1:ncycle){
      Mean <- Xj[[j]]%*%Mu[,k]
      denom <- (2*pi)^(-nchannel/2) * det( V[[k]] + Xj[[j]] %*% Sigma[[k]] %*% t(Xj[[j]]) )^(-1/2)
      SigmaInv <- solve( V[[k]] + Xj[[j]] %*% Sigma[[k]] %*% t(Xj[[j]]) )
      for(i in 1:nread){
        Z.denom[i,k,j] <- denom * exp(-0.5*t(Int[i,,j]-Mean) %*% SigmaInv %*% (Int[i,,j]-Mean) )
      }
    }
  }
  Z <- array(0,dim=c(nread,nchannel,ncycle))
  for(i in 1:nread){
    for(j in 1:ncycle){
      Z.sum <- sum(Z.denom[i,,j]*Pi)
      Z[i,,j] <- ( Z.denom[i,,j] * Pi ) / Z.sum
    }
  }

  a <- list(posterior.Z=Z,Pi=Pi,Mu=Mu,Sigma=Sigma,V=V,iter=iter,DensityInd=DensityInd)

  return(a)

} ## ends function
