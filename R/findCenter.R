#' A MVDA Function
#'
#' This function select cluster prototypes. It select the prototype as the more correlated with the other elements in the cluster
#' @param DB your matrix dataset
#' @param cluster.vector a numeric vector of clustering results
#' @keywords multi-view clustering; prototype; correlation
#' @return a matrix of prototypes
#' @export

findCenter <- function(DB,clust_vector){
    center <- NULL #matrice dei centroidi
    name_center <- c()
    elem_par_cluster <- table(clust_vector)
    for(i in names(elem_par_cluster)){
        # cat(i,"--> \n")
        if(elem_par_cluster[i]>1){
            clustElem <- DB[which(clust_vector==i),]
            matCor <- cor(t(clustElem), method="pearson");
            bestCenter <- which.max(apply(matCor,1,FUN=function(riga){
                sum(riga)
            }))
            center <- rbind(center,clustElem[bestCenter,])
            name_center <- c(name_center,attr(bestCenter,which="names"))
        }else{
            clustElem <- DB[which(clust_vector==i),]
            center <- rbind(center,clustElem)
            name_center <- c(name_center,rownames(DB)[which(clust_vector==i)])
        }
    }
    rownames(center)<-name_center
    return(center)
}