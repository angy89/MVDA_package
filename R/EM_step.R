EM_step <- function(vector.clust, DB,center ){
    new.vector.clust <- vector.clust;
    singleCluster <- c(which(table(vector.clust)==1),which(table(vector.clust)==2))
    notclustered_index <- c();
    for(i in 1:length(singleCluster)){
      notclustered_index <- c(notclustered_index,which(vector.clust==singleCluster[i]))
    }
    notClustered <-as.matrix(DB[notclustered_index,]); 
    for(i in 1:length(notclustered_index)){
      correlation <- apply(center[-singleCluster,],1,FUN=function(elem){
        cor(elem,notClustered[i,])
      })
      bestCenter <- which.max(correlation);
      new.vector.clust[notclustered_index[i]]<-bestCenter
    }
    tab<-table(new.vector.clust)
    for(i in 1:length(tab)){
      new.vector.clust[new.vector.clust == attr(tab[i],which="name")]<-i
    }
    EMcenter <- as.matrix(findCenter(DB,new.vector.clust));
    return(list(cluster = new.vector.clust,center = EMcenter))
}