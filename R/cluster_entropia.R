#' A MVDA Function
#'
#' This function evaluate the entropy of a cluster
#' @param cluster.vector a numeric vector of clustering results
#' @keywords multi-view clustering; entropy
#' @return the entropy in the cluster
#' @export
cluster_entropia<- function(vector.clust){
    tab <- table(vector.clust);
    tab <- tab/sum(tab)
    entropia <- sum(tab * log(tab+0.0001))
    return(-entropia)
}