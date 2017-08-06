#' A MVDA Function
#'
#' This function calculates the approximation error at each iteration as the difference between the matrix X and its approximation PH
#' @param X is the matrix resulting from the concatenation of the membership matrix for each view
#' @param PH is the approximation of X
#' @param l is the number of rows of X
#' @param n is the number of columns of X
#' @return the approximation error
#' @export

error <- function(X,PH,l,n){
  err <- 0

  for(i in 1:l){
    for(j in 1:n){
      err <- err + (X[i,j] - PH[i,j])^2
    }
  }

  return(err)
}

#' A MVDA Function
#'
#' This function calculates the matrix T indicating the contributing ot each view to each meta-cluster
#' @param P is the projection matrix of the approximation PH
#' @param k is the number of cluster
#' @param nView is the number of views in the experiments
#' @param V1.lc1 is the number of cluster in the first view
#' @param V2.lc2 is the number of cluster in the second view
#' @param V3.lc3 is the number of cluster in the third view and so on...
#' @return the approximation error
#' @export
find_T <- function(P,k,nView, V1.lc1,V2.lc2,V3.lc3,V4.lc4,V5.lc5){
  T <- matrix(0,nView,k);
  lc1 <- V1.lc1
  lc2 <- V2.lc2
  lc3 <- V3.lc3
  lc4 <- V4.lc4
  lc5 <- V5.lc5
  for(h in 1:nView){
    for(f in 1:k){
      if(h==1) T[h,f] <- (sum(P[1:lc1,f]))/sum(P[,f])
      if(h==2) T[h,f] <- (sum(P[(lc1+1):(lc1+lc2),f]))/sum(P[,f])
      if(h==3) T[h,f] <- (sum(P[(lc1+lc2+1):(lc1+lc2+lc3),f]))/sum(P[,f])
      if(h==4) T[h,f] <- (sum(P[(lc1+lc2+lc3+1):(lc1+lc2+lc3+lc4),f]))/sum(P[,f])
      if(h==5) T[h,f] <- (sum(P[(lc1+lc2+lc3+lc4+1):(lc1+lc2+lc3+lc4+lc5),f]))/sum(P[,f])

    }
  }
  return(T)
}

#' A MVDA Function
#'
#' This function calculates the best value of k through a series of executions (from kmin to kmax) of the algorithm of factorization in each step calculates the entropy of the matrix P of Responsibilities and chooses the partition with greater entropy value
#' @param kmin is the value of k from wich it start to search
#' @param kmax is the maximum value of k
#' @param X is the matrix resulting from the concatenation of the membership matrix for each view
#' @param eps is the precision parameter
#' @return the best factorization and the optimum value of k
#' @export

find_k <- function(kmin,kmax,X,eps){

  S_k <- NULL;
  Res_list <- list();
  for(index in kmin:kmax){
    cat(index,"\n")
    fattorizzazione_ottima <- matrix_factorization(X,index,eps,iter_max = 100);
    P <- fattorizzazione_ottima$P
    k_ottimo <- index
    P_norm <- apply(P,1,FUN=function(p){
      p <- p/sum(p)
    })

    P_norm[is.nan(P_norm)]<-0

    #Calcolo l'entropia
    entr <- apply(X=P_norm,1,FUN=function(p){
      s <- 0
      for(j in 1: length(p)){
        if(p[j]==0){
          s <- s+ p[j]
        }else{
          s <- s+ (p[j] * log(p[j]))
        }
      }
      (-1)/log(index) * s
    })

    s_k_ottimo <- 1-mean(unlist(entr))

    for(index2 in 1:100){
      cat(index2,"\n")
      pb <- txtProgressBar(min = 1, max = 100, initial= 1, style = 3);
      Res <- matrix_factorization(X,index,eps,iter_max = 10)
      P <- Res$P

      #Normalizzo P
      P_norm <- apply(P,1,FUN=function(p){
        p <- p/sum(p)
      })

      P_norm[is.nan(P_norm)]<-0

      entr <- apply(X=P_norm,1,FUN=function(p){
        s <- 0
        for(j in 1: length(p)){
          if(p[j]==0){
            s <- s+ p[j]
          }else{
            s <- s+ (p[j] * log(p[j]))
          }
        }
        (-1)/log(index) * s
      })

      s_k <- 1-mean(unlist(entr))
      if(s_k>s_k_ottimo){
        s_k_ottimo<-s_k;
        k_ottimo <- index
        fattorizzazione_ottima <- Res;
      }
      setTxtProgressBar(pb,index2);
    }
    S_k <- c(S_k,s_k_ottimo);
    Res_list <- lappend(Res_list,fattorizzazione_ottima);
    close(pb)

  }

  k <- which.max(S_k);
  cat("length Res_list ",length(Res_list),"\n")
  cat("length S_k: ",length(S_k),"\n")
  cat("k vale : ",k,"\n")
  return(list(factorization = Res_list[k][[1]], k = k_ottimo))

}

#' A MVDA Function
#'
#' This function implements the matrix factorization algorithm
#' @param X is the matrix resulting from the concatenation of the membership matrix for each view
#' @param k1 the desidered number of cluster
#' @param iter_mat the maximum allowed iteration
#' @param eps is the precision parameter
#' @return a list containing the two matri P and H
#' @export
#'

matrix_factorization <- function(X,k1,eps,iter_max){

  l <- dim(X)[1];
  n <- dim(X)[2];
  options(warn=-1)
    PH <- .nndsvd.internal(X, k1, flag=0)
  options(warn=0)
    
    P <- PH$W
    H <- PH$H

  PH <- P%*%H
  old_err = Inf
  err <- sqrt(error(X,PH,l,n)/(l*n))
  n_iter <-0;

  while(abs(old_err-err) > eps && n_iter <= iter_max){
    XHt <- X %*% t(H)
    PHHt <- P %*% H %*% t(H)
    P <- P * XHt / PHHt
    P[is.nan(P)]<-0;
    PtX <- t(P) %*% X
    PtPH <- t(P) %*% P %*% H
    H <- H * PtX / PtPH
    H[is.nan(H)]<-0
    PH <- P %*% H
    PH[is.nan(PH)]<-0;
    old_err = err;
    err <- sqrt(error(X,PH,l,n)/(l*n))
    n_iter <- n_iter +1;
  }

  toRet <- list(P = P, H=H)
}
