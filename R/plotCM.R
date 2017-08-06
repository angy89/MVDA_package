#' A MVDA Function
#'
#' This function plot a confusion matrix
#' @param CM is the confusion matrix
#' @param toSave if you want to save or not the image in a file. Default is FALSE
#' @param cartela is the path where the file will be saved
#' @keywords multi-view clustering; confusion matrix
#' @return a ggplot object
#' @export
#' 

plotCM <- function(CM,xlab="Clusters",ylab="Pam50",title="Patients Classes Confusion Matrix", toSave=F,cartella="./MV_clust/Kmeans/"){
  toPlot <- as.table(CM)
  require(reshape2)
  require(ggplot2)
  
  ggplot(melt(toPlot,c(xlab,ylab)), aes(Clusters,Pam50)) ->ggp
  ggp + geom_tile(aes(fill = value),colour = "white") ->ggp
  ggp + scale_fill_gradient(low = "aliceblue", high = "steelblue")->ggp
  ggp + ggtitle(title) -> ggp
  ggp
  
  if(toSave){
    cat(paste(cartella,gsub(" ","_", title),".jpg",sep=""))
    ggsave(filename=paste(cartella,gsub(" ","_", title),".jpg",sep=""),plot=ggp)
  }
  return(ggp)
}

