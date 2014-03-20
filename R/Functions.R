######################################################
# Hglasso Function for R package
######################################################
hgl <- function(S,lambda1,lambda2=100000,lambda3=100000,convergence=1e-10,maxiter=1000,start="cold",var.init=NULL,trace=FALSE){
	
  p <- nrow(S)
  rho <- 2.5
# Variables initialization
  if(start=="cold"){
  oldTheta <- Theta <- V <- Z <- diag(rep(1,p))
  W1 <- W2 <- W3 <- Gamma <- matrix(0,p,p)
  tildeTheta <- tildeV <- tildeZ <- matrix(0,p,p)
  }
  
  else if(start=="warm"){
  Theta <- diag(rep(1,p))	
  oldTheta <- var.init$Theta
  V <- var.init$V
  Z <- var.init$Z  
  W1 <- var.init$W1 
  W2 <- var.init$W2
  W3 <- var.init$W3
  Gamma <- var.init$Gamma  
  tildeTheta <- var.init$tildeTheta    
  tildeV <- var.init$tildeV
  tildeZ <- var.init$tildeZ
  }
  
  
  criteria <- 1e10 	
  i <- 1  	
# While loop for the iterations
  while(criteria > convergence && i <= maxiter){
  	
	if(trace==TRUE && i%%10==0) print(paste("iteration = ",i,sep=""))
	  
	Theta <- updateTheta(tildeTheta,W1,S,rho)
  	
  	Z <- updateZ(tildeZ,W3,lambda1,rho)	 

    V <- updateV(tildeV,W2,lambda2,lambda3,rho)

    Gamma <- updateGamma(Theta,V,Z,W1,W2,W3,rho)

    tildeTheta <- updatetildeTheta(Theta,W1,Gamma,rho)

    tildeV <- updatetildeV(V,W2,Gamma,rho)

    tildeZ <- updatetildeZ(Z,W3,Gamma,rho)
    
    W1 <- W1+Theta-tildeTheta
    
    W2 <- W2+V-tildeV
    
    W3 <- W3+Z-tildeZ
	
	Theta <- Z+V+t(V)	
	criteria <- sum((Theta-oldTheta)^2)/sum((oldTheta)^2)
	oldTheta <- Theta
	i <- i+1
	 	
  }  	
  Theta <- Z+V+t(V)

  if(i>maxiter){
  	 warning("The algorithm has not converged by the specified maximum number of iteration")
  }

  Obj <- HGLObj(Theta,S,Z,V,lambda1,lambda2,lambda3)
	
return(list(Theta=Theta,V=V,Z=Z,W1=W1,W2=W2,W3=W3,Gamma=Gamma,tildeTheta=tildeTheta,tildeV=tildeV,tildeZ=tildeZ,iteration=i,objective=Obj))	
}	

######################################################
# Objective Function for HGL
######################################################
HGLObj <- function(Theta,S,Z,V,lambda1,lambda2,lambda3){

tempV <- 0
p <- nrow(Theta)
A <- matrix(0,p,p)
B <- matrix(0,p,p)  
diag(A) <- diag(V)
diag(B) <- diag(Z)
tempV <- sum(sqrt(apply((V-A)^2,2,sum)))
  
return(-determinant(Theta,logarithm=TRUE)$modulus[1]+sum(S*Theta)+lambda1*sum(abs(Z-B))+lambda2*sum(abs(V-A))+lambda3*tempV)
}

######################################################
# Objective Function For Covariance Estimation
######################################################
HcovObj <- function(Sigma,S,Z,V,lambda1,lambda2,lambda3){

tempV <- 0
p <- nrow(Sigma)
A <- matrix(0,p,p)
B <- matrix(0,p,p)  
diag(A) <- diag(V)
diag(B) <- diag(Z)
tempV <- sum(sqrt(apply((V-A)^2,2,sum)))
  
return(0.5*sum((Sigma-S)^2)+lambda1*sum(abs(Z-B))+lambda2*sum(abs(V-A))+lambda3*tempV)
}

######################################################
# Soft-thresholding Operator
######################################################
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}

######################################################
# Update Theta
######################################################
updateTheta <- function(tildeTheta,W1,S,rho){
  C <- tildeTheta-W1-S/rho
  a <- eigen(C)
  D <- diag(a$values)
  U <- a$vectors
  Theta <- 1/2*U%*%(D+sqrt(D*D+4/rho*diag(rep(1,nrow(S)))))%*%t(U)		
  return(Theta)
}


######################################################
# Update Sigma
######################################################
updateSigma<-function(tildeSigma,W1,S,rho,epsilon){

temp1<-(1/(1+rho))*(S+rho*tildeSigma-rho*W1)
temp2<-eigen(temp1)
eigenvalue<-temp2$values
eigenvalue[which(eigenvalue<=1e-6)]<-epsilon
eigenvectors<-temp2$vectors
return(eigenvectors%*%diag(eigenvalue)%*%t(eigenvectors))
}

######################################################
# Update Z
######################################################
updateZ <- function(tildeZ,W3,lambda1,rho){
  A <- tildeZ-W3
  B <- lambda1/rho
  Z <- Soft(A,B)
  diag(Z) <- diag(A)
  return(Z)	
}

######################################################
# Update V
######################################################
updateV <- function(tildeV,W2,lambda2,lambda3,rho){
  p <- nrow(tildeV)
  V <- matrix(0,nrow=p,ncol=p)
  C <- tildeV-W2
  D <- C
  diag(D) <- 0
  for(j in 1:p){
    temp <- Soft(D[,j],lambda2/(rho))
    if(sum(temp^2)==0){
      V[,j] <- 0	
    }
    else if(sum(temp^2)!=0){
      V[,j] <- max(0,1-lambda3/(rho*(sqrt(sum(temp^2)))))*temp
    }
  }
  diag(V) <- diag(C)
  return(V)
}

######################################################
# Update tildeZ
######################################################
updatetildeZ <- function(Z,W3,Gamma,rho){
  tildeZ <- NULL
  tildeZ <- Gamma/rho+W3+Z
  return(tildeZ)	
}

######################################################
# Update tildeV
######################################################
updatetildeV <- function(V,W2,Gamma,rho){
  tildeV <- (Gamma+t(Gamma))/rho+V+W2
  return(tildeV)	
}

######################################################
# Update tildeTheta
######################################################
updatetildeTheta <- function(Theta,W1,Gamma,rho){
  tildeTheta <- W1+Theta-Gamma/rho
  return(tildeTheta)	
}

######################################################
# Update Gamma
######################################################
updateGamma <- function(Theta,V,Z,W1,W2,W3,rho){
  Gamma <- rho/6*(Theta+W1-V-t(V)-Z-W2-t(W2)-W3)
  return(Gamma)	
}
