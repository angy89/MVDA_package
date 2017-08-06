#' A MVDA Function
#'
#' This function concatenate list
#' @param lst list to be concatenated
#' @return a list
#' @export

lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
}