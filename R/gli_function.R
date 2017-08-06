#' A MVDA Function
#'
#' Function that perform che multiple view late integration clustering
#' @param A is a matrix in R^n*r (where n is the number of elements and r is the number of clusters) of m membership matrix stacked orizontally
#' @param k the desidered number of cluster 
#' @param alfa that is a positive constant
#' @param eps is the constant parameter for convergence
#' @return a list containing two matrix: - A clustering pattern matrix B (in R^ n*k) - A matrix P (in R^ k*r) of m mapping matrices stacked orizontally
#' @export
#' 


GLI <- function(A, k,alfa,eps){
  n <- dim(A)[1]
  r <- dim(A)[2]
  #Initilize B under the constarint B1 = 1 (the sum of all elements on each row of B must be 1)
  B <- matrix(data=runif(n=(n*k),min=1,max=10),nrow=n,ncol=k)
  B <- apply(B,MARGIN=1,FUN=function(row){
    return(row/sum(row))
  })
  B <- t(B)
  
  P <- matrix(data=runif(n=(k*r),min=1,max=10),nrow=k,ncol=r)
  gi_old <- GI(A,B%*%P)
  
  B <- Bupdate(A,B,P,alfa)
  P <- Pupdate(A,B,P)
  
  gi_new <- GI(A,B%*%P)
  
  niter <- 0
  while((gi_old  - gi_new) > eps){
    gi_old <- gi_new
    
    B <- Bupdate(A,B,P,alfa)
    P <- Pupdate(A,B,P)
    
    gi_new <- GI(A,B%*%P)
    niter <- niter + 1
  }
  cat("Number of iteration: ",niter,"\n")
  return(list(B = B, P=P))
}

#Update matrix P as step 4 of Algorithm 1 Multiple View Clustering
Pupdate <- function(A,B,Po){
  k <- dim(Po)[1]
  r <- dim(Po)[2]
  n <- dim(A)[1]
  
  P <- matrix(data=NA,nrow=k,ncol=r)
  
  for(g in 1:k){
    for(j in 1:r){
      
      num <- 0
      denum <- 0
      for(i in 1:n){
        num <- num + (A[i,j] * B[i,g] / sum((B[i,] * Po[,j])) )
        denum <- denum + B[i,g]
      }
      
      #Update
      P[g,j] <- Po[g,j] * (num / denum)
    }
  }
  return(P)
}

#Update matrix B as step 3 of Algorithm 1 Multiple View Clustering
Bupdate <- function(A,Bo,P,alfa){
  n <- dim(Bo)[1]
  k <- dim(Bo)[2]
  r <- dim(P)[2]
  B <- matrix(data=NA,nrow=n,ncol=k)
  for(i in 1:n){
    for(g in 1:k){
      #Find numnerator value
      s_num <- 0
      denum <- 0
      for(j in 1:r){
        s_num <- s_num + ( A[i,j] * (P[g,j]/ sum((Bo[i,] * P[,j])) ) )
        #Find Denumerator value
        denum <- denum + (P[g,j] + alfa)
      }
      
      num <- s_num + alfa / sum(Bo[i,])
      
      #Update value
      B[i,g] <- Bo[i,g] * (num / denum)
    }
  }
  
  return(B)
}

# Function that perform the Generalized I-divergence between two matrix X and Y
# Input:
#   - Two matrix X and Y of the same dimension
# Output:
#   - The value of Generalized I-divergence gi between X and Y
GI <- function(X,Y){
  X <- X + 0.1
  Y <- Y + 0.1
  nrow <- dim(X)[1]
  ncol <- dim(X)[2]
  
  gi <- 0
  
  for(i in 1:nrow){
    for(j in 1:ncol){
        gi <- gi + (log(X[i,j]) * log(X[i,j]/Y[i,j]) - X[i,j] + Y[i,j])
    }
  }
  
  return(gi)
}