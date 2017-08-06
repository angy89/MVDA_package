#' A MVDA Function
#'
#' This function evaluate clustering stability by performing multi-view integration perturbing the dataset via leave-one-out and comparing
#' NMI between different clustering results.
#' 
#' @param method is the method you can use for the integration. Values can be matrix_factorization (default) of gli.
#' @param A is the staked matrix of membeship matrix of each view.
#' @return a list containing the list of clustering obtained, the NxN matrix of NMI betwee each couple of clustering, the mean and the standard deviation of this matrix.
#' @export

mvclst <- function(method="MF",A,nClusters,nPatients, eps = 1e-6,iter_max = 100,alfa = 1 ){
  A <- as.matrix(A)
  # nClusters <- dim(A)[2]
  # nPatients <- dim(A)[1]
  pb <- txtProgressBar(min=1,max=nPatients,style=3)
  cluster_list <- list()
  indices = 1:nPatients
  
  if(method != "MF" & method!="GLI"){
    cat("Please specify one between the MF and GLI methods!")
    return(-1)
  }
  
  if(method =="MF"){
    A = t(A)
  }
  
  
  for(i in 1:nPatients){
    if(method =="MF"){
      matrix_factorization(X=A[,-i],k1=nClusters,eps=eps,iter_max=iter_max)->mf
      clustVector <- apply(mf$H,2,which.max)
    }
    if(method == "GLI"){
      GLI(A=A[-i,],k=nClusters,alfa=alfa,eps=eps)->GLI.res
      clustVector <- unlist(apply(GLI.res$B,1,which.max))
    }
    clustering = rep(NA,nPatients)
    clustering[indices[-i]] = clustVector
    cluster_list[[i]] = clustering
    setTxtProgressBar(pb=pb,value=i)
  }
  close(pb)  
  
  nmi_matrix <- matrix(0, nPatients, nPatients)
  pb = txtProgressBar(min = 1,max = nPatients,style=3)
  for(i in 1:nPatients){
    for(j in 1:nPatients){
      clus1 <- cluster_list[[i]]
      clus2 <- cluster_list[[j]]
      nas <- union(which(is.na(clus1)), which(is.na(clus2)))
      mutinfo <- mutinformation(clus1[-nas], clus2[-nas], method='emp')
      entr1 <- infotheo::entropy(clus1[-nas], method='emp')
      entr2 <- infotheo::entropy(clus2[-nas], method='emp')
      nmi_matrix[i,j] <- mutinfo / mean(entr1, entr2)
    }
    setTxtProgressBar(pb,i)
  } 
  close(pb)
  return(list(cluster_list = cluster_list,nmi_matrix = nmi_matrix, sd=sd(nmi_matrix),mean=mean(nmi_matrix)));  
}#end function

    
    