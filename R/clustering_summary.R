#' A MVDA Function
#'
#' This function evaluate the summary intdex for a clustering results
#' @param DB your matrix dataset
#' @param cluster is a numeric vector of clustering results
#' @keywords clustering evaluation
#' @return a vector of evaluation indexes
#' @export


clustering_summary <- function(DB, cluster){
  diss.cor <-1- cor(t(DB),method="pearson")
  toPrint <- matrix(0,nrow=1,ncol=6)
  NOBJS = dim(DB)[1]

  colnames(toPrint)<- c("#clusters","#singleton","Intra_avg","inter_complete","Compression_gain","Index")

    tab <- table(cluster)
    center = findCenter(DB = DB,clust_vector = cluster)
    j = 1
    toPrint[j,1]<-length(tab)

    toPrint[j,2]<-length(which(tab==1))


    cls.scatt.diss.mx(diss.mx=diss.cor,clust=cluster) -> clust_eval
    clust_eval$intracls.average -> intra.avg
    clust_eval$intercls.complete -> inter.c

    toPrint[j,3] = round(1-mean(intra.avg),digits = 2)
    toPrint[j,4] = round(1-mean(inter.c),digits = 2)
    toPrint[j,5] = round(1 - (length(table(cluster))/NOBJS),2)

    toPrint[j,6] =(1/3) * ( ((toPrint[j,3]+1)/2) + (1-(toPrint[j,4]+1)/2) +
                              (1-(toPrint[j,2]/length(table(cluster)))) ) # + toPrint[j,5])

    return(toPrint)
}
