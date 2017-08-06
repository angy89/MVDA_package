#' A MVDA Function
#'
#' This function plots the variance between of each features between the centroids
#' @param center is the matrix of centroids of the multi-view clusters; each row is a centroid andh the columns are the features
#' @return the correlation of the variables sorted in decreasing order
#' @export

plot_feature_variance_between_centroids = function(center,v_limit = 12,h_limit=1){
  var_col <- apply(mf_res$center,MARGIN=2,FUN=sd)
  var_col <- sort(var_col,decreasing=T)
  plot(var_col,type="l",main="Variance between cluster centroids",ylab="variance",xlab="",axes=FALSE)
  abline(h=h_limit,lty=2)
  abline(v=v_limit,lty=2)
  labs  <- names(var_col)
  text(x=1:length(labs),par("usr")[3], labels = labs, srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex=.5)
  axis(2)
  return(var_col)
}

#' A MVDA Function
#'
#' This function plots the correlation between the centroids
#' @param centers is the matrix of centroids of the multi-view clusters; each row is a centroid andh the columns are the features
#' @export
center_correlation = function(centers){
  corrgram(t(centers),order=TRUE, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt)
  
  # n_centers  = nrow(centers)
  # mat_cor <- matrix(0,n_centers,n_centers)
  # 
  # for(i in 1:n_centers){
  #   for(j in 1:n_centers){
  #     mat_cor[i,j] <- cor(centers[i,],centers[j,],method="pearson")        
  #   }
  # }
  # hist(mat_cor, main="correlazione tra i cluster")
}

#' A MVDA Function
#'
#' This function plots the boxplot of the features in the multi-view clusters. Each boxplot is coloured by the class labels assigned by the fisher test
#' @param centers is the matrix of centroids of the multi-view clusters; each row is a centroid andh the columns are the features
#' @param clusterMV is the vector of the multi view clustering
#' @param var_col are the feature ordered by variance across the centroids
#' @param nFeat is the number of feature to be considered for the boxplot
#' @param patientsDB is the dataset of the patient on the columns and the prototypes on the rows
#' @param cluster_class is the vector of the patient's classes assigned to each multi-view cluster by the fisher test
#' @param patClasses is the vector of the patients's classes
#' @param patClassColors is the vector of the color assigned to each class
#' @export
mv_clust_boxplot = function(centers, clusterMV,var_col,nFeat = 12,patientsDB,cluster_class,
                            patClasses=c("Basal","Her2","LumA","LumB","Normal"),
                            patClassColors = c("purple","blue","green","red","orange")){
  patientsDB= t(patientsDB)
  
  classi <- patClasses
  par(mfrow=c(4,3))
  n_centers = as.numeric(names(table(clusterMV)))
  
  if(nFeat>ncol(centers)) nFeat = ncol(centers)
  
  for(i in n_centers){
    if(sum(clusterMV==i)==1){
      patients.cl <- patientsDB[clusterMV==i,names(var_col)][1:nFeat]
      
    }else{
      patients.cl <- patientsDB[clusterMV==i,names(var_col)][,1:nFeat]
      
    }
    boxplot(patients.cl,xlab="",axes=FALSE,type="l",col = patClassColors[cluster_class[i]],
            ylab="Paziente",main=classi[cluster_class[i]])
    axis(2)
    labs  <- names(var_col)[1:12]
    text(x=1:12,par("usr")[3], 
         labels = labs, srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex=.5)
  }
  
}

#' A MVDA Function
#' 
#' This function plots centroids of each multi-view cluster
#' @param var_col are the feature ordered by variance across the centroids
#' @param nFeat is the number of features in the dataset
#' @param onlyVariables is a boolean if TRUE, only the most nvarFeat variables features will be used for the plot
#' @param nvarFeat is the number of features to be considered for the boxplot
#' @param cluster_class is the vector of the patient's classes assigned to each multi-view cluster by the fisher test
#' @param patClasses is the vector of the patients's classes
#' @export


# plot_centroids = function(cluster_class,patClasses,mv_clust,
#                           nFeat = 12,centers,var_col,nvarFeat = 6,onlyVariables=T){
#   cluster_class=cluster_class[as.numeric(names(table(mv_clust)))]
#   if(onlyVariables){
#     centers = centers[,names(var_col)[1:nvarFeat]]
#     Labels<-patClasses[cluster_class]
#     #Labels <- Labels[-9]
#     Labels<-sapply(Labels, function (x) rep(x,nvarFeat))
#     Labels<-as.vector(Labels)
#   }else{
#     Labels<-patClasses[cluster_class]
#     #Labels <- Labels[-9]
#     Labels<-sapply(Labels, function (x) rep(x,nFeat))
#     Labels<-as.vector(Labels)
#   }
#   
#   
#   sa <- stack(as.data.frame(t(centers)))
#   sa$x <- rep(seq_len(ncol(centers)), nrow(centers))
#   qplot(x, values, data = sa, colour = Labels, geom = "line") +
#     xlab(" ")+
#     ylab(" ")+
#     labs(colour = "Centroids") +
#     scale_x_continuous(breaks=1:12,labels=colnames(cluster_class))+
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   sa$label <- Labels
#   ggplot(sa, aes(x, values, colour=as.factor(label))) + 
#     geom_line() + 
#     facet_wrap(~label,nrow=5) +
#     xlab(" ")+
#     ylab(" ")+
#     labs(colour = "Centroids") +
#     scale_x_continuous(breaks=1:12, 
#                        labels=colnames(centers))+
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   
# }
# 
