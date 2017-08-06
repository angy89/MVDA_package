#' A MVDA Function
#'
#' This function execute partitioning around medoids algoritm using function from fpc package.
#' 
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @param initial_medoids index of the desidered initial medoids you want to use to initialize the algorithm
#' @keywords pamk; ward
#' @return a list containing three field: pamk.res is the pamk results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export
#' 

pamk_preprocessing<-function(DB=NULL,nCenters = NULL,initial_medoids=NULL){
  require(fpc)
  if(is.null(DB)){
    stop("You must insert a DB!\n")
  }else{
    if(is.null(nCenters)){
      stop("You must insert the number of centers!\n")
    }
    else{
      cat("Execute pamk clustering...\n")
      DB.diss <- as.dist(1-cor(x=t(DB),method="pearson",use="pairwise.complete.obs"));
      
      pamk(data=DB.diss,krange = nCenters,diss=TRUE,medoids=initial_medoids)->DB_pamk
      DB_pamk$pamobject$clustering -> clusters
      DB[DB_pamk$pamobject$medoids,] -> centers
      toRet <- list(pamk.result = DB_pamk, clustering = clusters, centers = centers);
      cat("End...\n")
      return(toRet);
    }
  }
}
