#' A MVDA Function
#'
#' This function find correlation between object in each cluster and its prototype
#' @param DB your matrix dataset
#' @param cluster.vector a numeric vector of clustering results
#' @param center a matrix with cluster prototype
#' @keywords multi-view clustering; feature-relevance; correlation
#' @return the vector of mean correlation between objects in the cluster and its prototypes
#' @export
object_correlation <- function(DB, clust_vector,center){
    myApply <- lapply
    DB <- as.matrix(DB)
    center <- as.matrix(center)
    corr_between_object_and_centroid <- c();
    
    for(i in 1:dim(DB)[1]){
        corr_between_object_and_centroid[i] <- cor(DB[i,],center[clust_vector$x[i],])
    }
    
    return(corr_between_object_and_centroid)
}