#' A MVDA Function
#'
#' This function print statistics about clustering results
#' @param DB your matrix dataset
#' @param cluster.vector a numeric vector of clustering results
#' @param center a matrix with cluster prototype
#' @param folder is the path where the function will save the statistics files
#' @keywords multi-view clustering; feature-relevance; correlation
#' @export
printStatistics<-function(DB, clust_vector,center,folder){
    library("gdata")
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
    toPrint <- cbind(corr_between_object_and_centroid,table(clust_vector))
    toPrint <- toPrint[order(toPrint[,1]),]
    rownames(toPrint) <- names(table(clust_vector))
    colnames(toPrint) <- c("Correlation","ElemInCluster")
    write.table(toPrint,paste(folder,"correlazioni_centroidi_cluster.txt",sep="\t"))
    pdf(file=paste(folder,"/MeanCorrelation.pdf",sep=""))
    hist(corr_between_object_and_centroid,breaks=60,freq=F,main="Mean correlation between prototypes and elements in cluster",xlab="Correlazione")
    
    lines(density(corr_between_object_and_centroid),col="black")
    dev.off();
    cat("Mean correlation: ",mean(corr_between_object_and_centroid),"\n");
    
    inCorrelation <- list()
    for(i in 1:n_cluster){
        elem.i <- which(clust_vector==attr(elem_par_cluster[i],which="name"))
        if(length(elem.i)==1){
            inCorrelation <- lappend(inCorrelation,1)
        }else{
            inCorrelation <- lappend(inCorrelation,cor(t(DB[elem.i,])))
        }
    }
    mean_correlation <- myApply(inCorrelation, FUN=function(mi){
        if(class(mi)=="numeric")return(1)
        mean(lowerTriangle(mi, diag=FALSE))
    })
    
    mean_correlation <- unlist(mean_correlation);
    #dev.new()
    #X11()
    pdf(file=paste(folder,"/AverageCorrelation.pdf",sep=""))
    
    hist(mean_correlation,freq=F,main="Mean correlation in clusters",xlab="Correlation")
    lines(density(mean_correlation),col="black")
    dev.off()
    sd_correlation <- myApply(inCorrelation, FUN=function(mi){
        #sqrt(var(lowerTriangle(mi, diag=FALSE)))
        var(lowerTriangle(mi, diag=FALSE))
    })  
    sd_correlation <- unlist(sd_correlation);
    na <- which(is.na(sd_correlation))
    sd_correlation[na]<-0
    #dev.new()
    #x11()
    pdf(file=paste(folder,"/SDCorrelation.pdf",sep=""))
    
    hist(sd_correlation,freq=F, breaks=31, main="Mean standard deviation in cluster",xlab="Standard deviation")
    lines(density(sd_correlation),col="black")
    dev.off()
}