#' A MVDA Function
#'
#' This function find correlation between object in each cluster and its prototype
#' @param DB your matrix dataset
#' @param cluster.vector a numeric vector of clustering results
#' @param center a matrix with cluster prototype
#' @keywords multi-view clustering; feature-relevance; correlation
#' @return the vector of mean correlation between object and prototypes
#' @export
mean_correlation <- function(DB, clust_vector,center){
    myApply <- lapply
    DB <- as.matrix(DB)
    center <- as.matrix(center)
    corr_between_object_and_centroid <- c();
    
    elem_par_cluster <- table(clust_vector)
    n_cluster <- length(elem_par_cluster)
    for(i in 1:n_cluster){ #per ogni cluster
        elemClust.i <- which(clust_vector==attr(elem_par_cluster[i],which="name"))
        mean_corr <- c()
        for(j in 1:length(elemClust.i)){
            mean_corr <- c(mean_corr,cor(center[i,],DB[elemClust.i[j],],method="pearson"))  
        }
        mean_corr <- mean(mean_corr)
        corr_between_object_and_centroid <-c(corr_between_object_and_centroid,mean_corr);
    }
    
    return(corr_between_object_and_centroid)
}