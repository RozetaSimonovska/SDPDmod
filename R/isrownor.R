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
#' @examples
#' data("usa46", package="SDPDmod")
#' isrownor(usa46)
#'
#' @seealso \code{\link{rownor}}
#
#' @export

isrownor<-function(W){

  if(nrow(W)!=ncol(W)) stop("Error in matrix!")

  if(isTRUE(all.equal(as.vector(rowSums(W)), 
                      rep(1, nrow(W)),
                      tol = .Machine$double.eps^0.5))) {
    return(TRUE)} else {
      return(FALSE)}

}
