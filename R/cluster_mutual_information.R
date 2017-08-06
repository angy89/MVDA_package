#' A MVDA Function
#'
#' This function evaluate the mutual information between two clustering results
#' @param vc1 a numeric vector of clustering results
#' @param vc2 a numeric vector of clustering results
#' @param elem the matrix representing your dataset
#' @keywords multi-view clustering; entropy
#' @return the mutual information between clustering 
#' @export

cluster_mutual_information<- function(vc1,vc2,elem){
    tab1 <- table(vc1);
    tab2 <- table(vc2);
    N <- sum(tab1)
    c <- 0;
    
    for(k in 1:length(tab1)){
        for(j in 1:length(tab2)){
            elem.ck <- elem[which(vc1==k)];
            elem.cj <- elem[which(vc2==j)];
            sum(elem.ck %in% elem.cj) -> nelem_intersect
            wk <- tab1[k]
            cj <- tab2[j]
            if(nelem_intersect == 0 ){
                value <- 0  
            }else{
                value <- (nelem_intersect/N) * log((N*nelem_intersect)/(wk*cj))                
            } 
            c <- c + value
        }
    }
    return(c)
}