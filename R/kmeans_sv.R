#' A MVDA Function
#'
#' This function execute kmeans clustering on single view patient prorotypes. It require library amap.
#' 
#' @param nCenters is the number of cluster we want to obtain
#' @param prototype is the matrix of prototype we want to cluster
#' @param method is the method by wich distance is evaluated. Default is pearson.
#' @param corr.use is an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". Default value="pairwise.complete.obs" 
#' @param iter.max The maximum number of iterations allowed to Kmeans
#' @param nstart If nCenter is a number, how many random sets should be chosen?
#' @keywords kmeans-clustering; ward
#' @return a list containing three field: pamk.res is the pamk results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

kmeans_sv <- function(nCenters,prototype,method="pearson",
                    corr.use="pairwise.complete.obs",iter.max=10000,nstart=10){
  
  require(amap)
  if(is.null(prototype) || is.null(nCenters)){
    stop("You must insert prototypes and the number of desidered clusters!\n")
  }else{
    N_PAT <- dim(prototype)[2]
    cat("Calculating dissimilarity matrix...\n")
    
    Kmeans(t(prototype),centers = nCenters,iter.max=iter.max,nstart = nstart,method = method)->DB_pamk
    cluster.m <- DB_pamk$cluster
    
    center <- findCenter(DB=prototype,clust_vector=cluster.m)
    toRet <- list(pamk.res=DB_pamk,clustering = cluster.m, center= center)
    return(toRet)
  }
}