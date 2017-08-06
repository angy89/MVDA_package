#' A MVDA Function
#'
#' This function execute hierarchical clustering with ward method on single view patient prorotypes
#' 
#' @param nCenters is the number of cluster we want to obtain
#' @param prototype is the matrix of prototype we want to cluster
#' @param corr_method is the method by wich distance is evaluated. Default is pearson.
#' @param use is an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". Default value="pairwise.complete.obs" 
#' @param hls_method is the the agglomeration method to be used for the hclust function. Default value is "ward".
#' @keywords hierarchical-clustering; ward
#' @return a list containing three field: hls.res is the hclust results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export
ward_sv <- function(nCenters,prototype,corr_method="pearson",
                    use="pairwise.complete.obs",hls_method="ward"){
  
  if(is.null(prototype) || is.null(nCenters)){
    stop("You must insert prototypes and the number of desidered clusters!\n")
  }else{
      N_PAT <- dim(prototype)[2]
      cat("Calculating dissimilarity matrix...\n")
      DB.dist <- as.dist(1-cor(x=prototype,method=corr_method,use=use));
      DB.dist[is.na(DB.dist)] <- 0
      cat("Calculating hierarchical clustering...\n")
      hls <- hclust(d=DB.dist^2,method=hls_method)
      cluster.m <- cutree(hls,nCenters)
      center <- findCenter(DB=prototype,clust_vector=cluster.m)
      
      toRet <- list(hls.res=hls,clustering = cluster.m,center=center)
      return(toRet)
  }
}
