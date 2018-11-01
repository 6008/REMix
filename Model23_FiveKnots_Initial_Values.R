initial <- function(int,sampsize){
  readind <- sample(1:dim(int)[1],sampsize)
  newint <- int[readind,,]
  nread <-  dim(newint)[1]
  nchannel <- dim(newint)[2]
  ncycle <- dim(newint)[3]
  l <- 1
  first <- second <- third <- fourth <- 1:(nread*ncycle)
  for(r in 1:nread){
    for(j in 1:ncycle){
      first[l] <- sort(newint[r,,j])[4]
      second[l] <- sort(newint[r,,j])[3]
      third[l] <- sort(newint[r,,j])[2]
      fourth[l] <- sort(newint[r,,j])[1]
      l <- l + 1
    }
  }
  X <- rep(0:(ncycle-1),nread)
  D = 1
  K = 5
  knots = quantile(1:nread,c(0.1667,0.3333,0.50,0.6667,0.8333))
  X1 = outer(X,1:D,"^")
  X2 = outer(X,knots,">") * outer(X,knots,"-")^D
  myX <- cbind(X1,X2)
  first_reg <- lm(first~myX)
  second_reg <- lm(second~myX)
  third_reg <- lm(third~myX)
  fourth_reg <- lm(fourth~myX)
  Mu <- matrix(0,ncol=nchannel,nrow=28)
  Mu[,1] <- c(first_reg$coef[1],second_reg$coef[1],third_reg$coef[1],fourth_reg$coef[1],first_reg$coef[2],second_reg$coef[2],third_reg$coef[2],fourth_reg$coef[2],first_reg$coef[3],second_reg$coef[3],third_reg$coef[3],fourth_reg$coef[3],first_reg$coef[4],second_reg$coef[4],third_reg$coef[4],fourth_reg$coef[4],first_reg$coef[5],second_reg$coef[5],third_reg$coef[5],fourth_reg$coef[5],first_reg$coef[6],second_reg$coef[6],third_reg$coef[6],fourth_reg$coef[6],first_reg$coef[7],second_reg$coef[7],third_reg$coef[7],fourth_reg$coef[7])
  Mu[,2] <- c(second_reg$coef[1],first_reg$coef[1],fourth_reg$coef[1],third_reg$coef[1],second_reg$coef[2],first_reg$coef[2],fourth_reg$coef[2],third_reg$coef[2],second_reg$coef[3],first_reg$coef[3],fourth_reg$coef[3],third_reg$coef[3],second_reg$coef[4],first_reg$coef[4],fourth_reg$coef[4],third_reg$coef[4],second_reg$coef[5],first_reg$coef[5],fourth_reg$coef[5],third_reg$coef[5],second_reg$coef[6],first_reg$coef[6],fourth_reg$coef[6],third_reg$coef[6],second_reg$coef[7],first_reg$coef[7],fourth_reg$coef[7],third_reg$coef[7])
  Mu[,3] <- c(third_reg$coef[1],fourth_reg$coef[1],first_reg$coef[1],second_reg$coef[1],third_reg$coef[2],fourth_reg$coef[2],first_reg$coef[2],second_reg$coef[2],third_reg$coef[3],fourth_reg$coef[3],first_reg$coef[3],second_reg$coef[3],third_reg$coef[4],fourth_reg$coef[4],first_reg$coef[4],second_reg$coef[4],third_reg$coef[5],fourth_reg$coef[5],first_reg$coef[5],second_reg$coef[5],third_reg$coef[6],fourth_reg$coef[6],first_reg$coef[6],second_reg$coef[6],third_reg$coef[7],fourth_reg$coef[7],first_reg$coef[7],second_reg$coef[7])
  Mu[,4] <- c(fourth_reg$coef[1],third_reg$coef[1],second_reg$coef[1],first_reg$coef[1],fourth_reg$coef[2],third_reg$coef[2],second_reg$coef[2],first_reg$coef[2],fourth_reg$coef[3],third_reg$coef[3],second_reg$coef[3],first_reg$coef[3],fourth_reg$coef[4],third_reg$coef[4],second_reg$coef[4],first_reg$coef[4],fourth_reg$coef[5],third_reg$coef[5],second_reg$coef[5],first_reg$coef[5],fourth_reg$coef[6],third_reg$coef[6],second_reg$coef[6],first_reg$coef[6],fourth_reg$coef[7],third_reg$coef[7],second_reg$coef[7],first_reg$coef[7])
  Sigma <- list()
  Sigma[[1]] <- diag(28)
  Sigma[[2]] <- diag(28)
  Sigma[[3]] <- diag(28)
  Sigma[[4]] <- diag(28)
  Pi <- rep((1/nchannel),nchannel)
  V <- list()
  V[[1]] <- V[[2]] <- V[[3]] <- V[[4]] <- diag(4)
  iv <- list(Mu,Sigma,Pi,V)
  names(iv) <- c("Mu","Sigma","Pi","V")
  return(iv)
}