#' @name isrownor
#'
#' @title Is the matrix row-normalized
#'
#' @description Checks if a spatial weights matrix is row-normalized.
#'
#' @param W spatial weights matrix
#'
#' @return Logical value. If the weights matrix is row-normalized
#' such that all rows sum up to 1, the value is TRUE.
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{rownor}}
#
#' @export

isrownor<-function(W){

  if(nrow(W)!=ncol(W)) stop("Error in matrix!")

  if(all(rowSums(W)==1)) {return(TRUE)} else {return(FALSE)}

}
