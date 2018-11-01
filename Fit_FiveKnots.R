qtrim <- 0.01
qknot <- c(0.1667,0.3333,0.50,0.6667,0.8333)
knots <- quantile(0:101,qknot)
library(ShortRead)
obj_id <- paste("Model23_FiveKnots_",Name,"_trim",qtrim,sep="")
load("Int_HQ100.rda")
load("SrfimEstimates_HQ100.rda")
source("Model23_FiveKnots_Initial_Values.R")
nread <- dim(int)[1]
nchannel <- dim(int)[2]
ncycle <- dim(int)[3]
source("/bigdata/cuilab/shared/basecall/Model23_FiveKnots_Initial_Values.R")
init <- initial(int=int,sampsize=100) # calculating initial values
Mu <- init$Mu
Sigma <- init$Sigma
V <- init$V
Pi <- init$Pi
niter <- 10
source("MCEM_MultipleKnots.R")
ptm1 <- proc.time()
output <- MCEM(Int=int,Mu=Mu,Sigma=Sigma,V=V,Pi=Pi,verb=TRUE,niter=niter,probs=probs,knots=knots,q=qtrim)
ptm2 <- proc.time()
TotTime <- ptm2 - ptm1
save.image(paste(obj_id,"_EM.rda",sep=""))
ReadInd <- output$DensityInd[,1]
CycleInd <- output$DensityInd[,2]
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
Ker <- matrix(0,ncol=nchannel,nrow=length(ReadInd))
for(l in 1:length(ReadInd)){
  i <- ReadInd[l]
  j <- CycleInd[l]
  for(k in 1:nchannel){
    Ker[l,k] <- t(int[i,,j] - Xj[[j]]%*%output$Mu[,k]) %*% solve(output$V[[k]]+Xj[[j]]%*%output$Sigma[[k]]%*%t(Xj[[j]])) %*% (int[i,,j] - Xj[[j]]%*%output$Mu[,k])
   }
}
zseq <- matrix(NA,nrow=nread,ncol=ncycle)
post <- output$posterior.Z
for(l in 1:length(ReadInd)){
  post[ReadInd[l],,CycleInd[l]] <- -Ker[l,]
}
for(i in 1:nread){
  for(j in 1:ncycle){
    nc <- which(post[i,,j]==max(post[i,,j]))
    if(nc==1) zseq[i,j]<-"A"
    if(nc==2) zseq[i,j]<-"C"
    if(nc==3) zseq[i,j]<-"G"
    if(nc==4) zseq[i,j]<-"T"
  }
}
z_seq <- apply(zseq,1,function(x) paste(x,sep="",collapse=""))
samp <- readFastq(paste("/bigdata/cuilab/shared/basecall/phiX174_HiSeq_Datasets/",Name,".fastq",sep=""))
Model23 <- ShortReadQ(sread=DNAStringSet(z_seq), quality=samp@quality, id=samp@id)
writeFastq(Model23,compress=FALSE,paste(obj_id,".fastq",sep=""))
source("/bigdata/cuilab/shared/basecall/FastqToAln.R")
aln(fasta_file="/bigdata/cuilab/shared/phiX174/phiX174",fastq_file=obj_id)
source("/bigdata/cuilab/shared/basecall/mismatch_flags.R")
Results <- mismatches(fasta_file="/bigdata/cuilab/shared/phiX174/phiX174.fasta",bam_file=paste(obj_id,".bam",sep=""),obj_id=obj_id,tile="")
write.table(round(t(Results),4),file=paste(obj_id,"_Results.txt",sep=""),quote=FALSE,row.names=FALSE,sep=", ")
Results
