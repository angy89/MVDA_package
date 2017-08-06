#' A MVDA Function
#'
#' This function execute kmeans clustering on single view patient prorotypes. It require library amap.
#' 
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @param method is the method by wich distance is evaluated. Default is pearson.
#' @param iter.max The maximum number of iterations allowed to Kmeans
#' @param nstart If nCenter is a number, how many random sets should be chosen?
#' @keywords kmeans-clustering; ward
#' @return a list containing three field: pamk.res is the pamk results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

kmeans_preprocessing<-function(DB=NULL,nCenters = NULL,iter.max=100,nstart=10,method="pearson"){
  require(amap)
  if(is.null(DB)){
    stop("You must insert a DB!\n")
  }else{
    if(is.null(nCenters)){
      stop("You must insert the number of centers!\n")
    }
    else{
      cat("Execute kmeans clustering...\n")
      Kmeans(DB,centers = nCenters,iter.max=iter.max,nstart = nstart,method = method)->DB_km
      DB_km$cluster -> clusters
      cat("Find centers...\n")
      findCenter(DB,clusters)->centers
      toRet <- list(km.result = DB_km, clustering = clusters, centers = centers);
      cat("End...\n")
      return(toRet);
    }
  }
}

