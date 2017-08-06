#' A MVDA Function
#'
#' This function calculate clustering by using rectangular som in bootrapping mode. We construct a similarity matrix of the sample by using SOM clustering results that is the imput for hierarchical clustering. 
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @param nboot is the number of iteration for the bootstrap method
#' @param xdim is the number of unit on the x-axis
#' @param ydim is the number of unit on the y-axis
#' @param t is a threshold for the consensus matrix. value below t are set to 0, and value above t are set to 1.
#' @description Given the matrix of prototype we execute clustering with SOM in bootstrap mode. In each bootstrap iterarion we permuted the dataset and trained a rectangular SOM. We then count how many time each couple of object have been clustered together normalized for the number of time they have been selected together in the dataset. We construct a matrix with this information for each couple of object that is used as distance for a hierarchical clustering with ward method. 
#' @keywords kohonen-SOM; correlation; hierarchical-clustering; ward
#' @return a list containing three field: som.res is the hierarchical clustering results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

SOM_preprocessing<-function(DB=NULL,xdim=NULL,ydim=NULL,nboot=100, t = 0.95){
  if(is.null(DB)){
    stop("You must insert a DB!\n")
  }else{
    if(is.null(xdim) || is.null(ydim)){
      stop("You must insert xdim and ydim!\n")
    }
    else{
      cat("Execute clustering with SOM...\n")
      M <- matrix(0,nrow=nrow(DB),ncol=nrow(DB));
      I <- matrix(0,nrow=nrow(DB),ncol=nrow(DB)); #la prima volta sono stati scelti tutti gli oggetti quindi sono presenti tutte le possibili coppie (i,j)
      rownames(M)<- colnames(M)<-paste("m",1:nrow(DB),sep="");
      colnames(I)<-rownames(I)<-1:nrow(DB);
      cat("executing bootstrap...\n")
      pb <- txtProgressBar(min=1,max=nboot,style=3);
      for(j in 1:nboot){
        smpl <- sort(sample(1:nrow(DB),size=nrow(DB),replace=TRUE))
        som.koh.obj <- som(as.matrix(DB[smpl,]),grid=somgrid(xdim=xdim,ydim=ydim,topo="rectangular"))
        clust <- som.koh.obj$unit.classif;
        elem_par_cluster <- table(clust)
        nclust <- length(elem_par_cluster)
        for(l in 1:nclust){
          index <- smpl[which(clust==as.integer(attr(elem_par_cluster[l],which="name")))] #indici degli elementi che sono stati raggruppati nello stesso cluster
          M[index,index]<-M[index,index]+1
        }
        I[smpl,smpl] <- I[smpl,smpl] + 1;  
        gc();
        setTxtProgressBar(pb=pb,value=j)
      }
      close(pb);
      M2 <- M/I;
      M2[is.nan(M2)] <- 0;
      gc()
      M2<-ifelse(M2>=t,1,0);
      diag(M2)<-0
      M2diss <- 1-M2;
      gc();
      M2diss <- M2diss * t(M2diss)
      gc();
      hls <- hclust(as.dist(M2diss)^2,method="ward.D");
            
      hls.clust <- cutree(hls,k=(xdim * ydim));
      center.hls <- as.matrix(findCenter(DB,hls.clust));
      
      toRet <- list(hls.result = hls, clustering = hls.clust, centers = center.hls);
      cat("End...\n")
      return(toRet);
    }
  }
}

