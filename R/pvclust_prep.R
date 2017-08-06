#' A MVDA Function
#'
#' This function calculate clustering of single view patient prototypes by using pvclust packages.
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @param hclust.method is the the agglomeration method to be used for the hclust function. Default value is "ward".
#' @param nboot is the number of bootstrap replications. The default is 100.
#' @param r is numeric vector which specifies the relative sample sizes of bootstrap replications. The default is 1.
#' @param dist.method is the distance measure to be used.
#' @param nCenters is the desidered number of clusters
#' @keywords pvclust; correlation; hierarchical-clustering; ward; p-value
#' @return a list containing three field: pvclust.res is the pvclust clustering results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

pvclust_prep<-function(DB,hclust.method="ward",nboot=100,
                       dist.method="correlation",alpha=0.95,
                       show.print=F,r=1,printLocation=".",km_center = F,nCenters){
  require(pvclust)
  if(is.null(DB)){
    stop("You must insert a DB!\n")
  }else{
    if(is.null(nCenters)){
      stop("You must insert the number of centers!\n")
    }
    else{
      cat("Execute pvcluster...\n")
      pvclust.res <- pvcluster(DB,hclust.method,nboot,dist.method,alpha,show.print,r,printLocation)
      pvclust.res$n_cluster -> n_cluster
      pvclust.res$cluster.vector -> cluster.vector
      cat("Finding centers...\n")
      center <- findCenter(DB=t(DB),clust_vector=cluster.vector)
      
      tab.cluster <- table(cluster.vector)
      mean(tab.cluster) -> mean.cl
      max(tab.cluster) -> max.cl
      min(tab.cluster)->min.cl
      
      if(min.cl == 1 || min.cl == 2){
        res <- EM_step(vector.clust=cluster.vector,DB=t(DB),center = center)
        cluster.vector <- res$cluster
        center <- findCenter(DB = t(DB),clust_vector=cluster.vector)
      }
      
      toRet <- list(pvclust.res = pvclust.res, clustering = cluster.vector, centers = center,n_cluster = n_cluster);
      cat("End...\n")
      return(toRet);
    }
  }
}

