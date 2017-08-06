#' A MVDA Function
#'
#' This function calculate clustering of single view patient prototypes by using pvclust packages.
#' @param prototype is the matrix of prototype we want to cluster
#' @param hls_method is the the agglomeration method to be used for the hclust function. Default value is "ward".
#' @param nboot is the number of bootstrap replications. The default is 100.
#' @param r is numeric vector which specifies the relative sample sizes of bootstrap replications. The default is 1.
#' @param dist.method is the distance measure to be used.
#' @keywords pvclust; correlation; hierarchical-clustering; ward; p-value
#' @return a list containing three field: pvclust.res is the pvclust clustering results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

pvclust_sv <- function(prototype,hls_method="ward",nboot=100,
                       t=0.95,print=F,r=1,print_folder=".",dist.method="correlation"){
  require(pvclust)
  if(is.null(prototype)){
    stop("You must insert prototypes and the number of desidered clusters!\n")
  }else{
    N_PAT <- dim(prototype)[2]
    cat("Calculating dissimilarity matrix...\n")
    
    pvclust.res <- pvcluster(prototype,hls.method,nboot,dist.method,t,print,r,print_folder)
    pvclust.res$cluster.vector -> cluster.m 
    
    center <- findCenter(DB=prototype,clust_vector=cluster.m)
    
    toRet <- list(pvclust.res=pvclust.res,clustering = cluster.m,center=center)
    return(toRet)
  }
}
