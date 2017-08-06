#' A MVDA Function
#'
#' This function execute hierarchical clustering with ward method
#' 
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @keywords hierarchical-clustering; ward
#' @return a list containing three field: hls.res is the hclust results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export
#' 
ward_hls_prep<-function(DB=NULL,nCenters = NULL){
    if(is.null(DB)){
      stop("You must insert a DB!\n")
    }else{
      if(is.null(nCenters)){
        stop("You must insert the number of centers!\n")
      }
      else{
        cat("Execute hierarchical clustering with ward methods...\n")
        DB.dist <- as.dist(1-cor(x=t(DB),method="pearson",use="pairwise.complete.obs"));
        hls <- hclust(d=DB.dist^2,method="ward.D")
        hls.clust <- cutree(hls,nCenters)
        cat("Finding centers...\n")
        hls.center <- as.matrix(findCenter(DB,hls.clust));
        toRet <- list(hls.result = hls, clustering = hls.clust, centers = hls.center);
        cat("End...\n")
        return(toRet);
      }
    }
}
  