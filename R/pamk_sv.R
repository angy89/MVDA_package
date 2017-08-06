#' A MVDA Function
#'
#' This function execute partitioning around medoids clustering with on single view patient prorotypes. It require library fpc.
#' 
#' @param nCenters is the number of cluster we want to obtain
#' @param prototype is the matrix of prototype we want to cluster
#' @param dist.method is the method by wich distance is evaluated. Default is pearson.
#' @param corr.use is an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". Default value="pairwise.complete.obs" 
#' @keywords pamk-clustering; ward
#' @return a list containing three field: pamk.res is the pamk results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

pamk_sv <- function(nCenters,prototype,dist.method="pearson",corr.use="pairwise.complete.obs"){
  
  require(fpc)     
  if(is.null(prototype) || is.null(nCenters)){
    stop("You must insert prototypes and the number of desidered clusters!\n")
  }else{
    N_PAT <- dim(prototype)[2]
    cat("Calculating dissimilarity matrix...\n")
    
    DB.diss <- 1-cor(x=prototype,method=dist.method,use=corr.use);
    DB.diss[is.na(DB.diss)] <- 0
    cat("Executing pamk... \n")
    pamk(data=DB.diss,krange = nCenters,diss=TRUE)->DB_pamk
    cluster.m <- DB_pamk$pamobject$clustering
    
    center <- findCenter(DB=prototype,clust_vector=cluster.m)
    toRet <- list(pamk.res=DB_pamk,clustering = cluster.m, center= center)
    cat("end. \n")
    return(toRet)
  }
}
