#' @name isrownor
#'
#' @title Is the matrix row normalised
#'
#' @description Check if a spatial weights matrix is row normalised.
#'
#' @param W spatial weights matrix
#'
#' @return Logical value. If the weights matrix is row normalised
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
