#' A MVDA Function
#'
#' This function consensus matrix from a list of clustering results. It evaluetes for each couple of patients
#' how many times they have been clustered together.
#' @param cluster_list is list of vector of clustering results
#' @keywords multi-view clustering; confusion matrix
#' @return a consensus matrix
#' @export
#' 
consensus_matrix <- function(cluster_list){
    cm <- matrix(0,nrow=length(cluster_list[[1]]),ncol=length(cluster_list[[1]]))
    cl.matrix <- matrix(unlist(cluster_list),nrow=length(cluster_list),ncol=length(cluster_list[[1]]),byrow=TRUE)
    
    #Per ogni coppia di pazienti
    for(i in 1:dim(cm)[1]){
        for(j in 1:dim(cm)[2]){
            cm[i,j] <- sum(cl.matrix[,i] == cl.matrix[,j])
        }
    }
    
    cm <- cm / dim(cl.matrix)[1];
    return(cm)

}