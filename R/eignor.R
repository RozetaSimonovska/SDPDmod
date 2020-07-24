#' @name eignor
#' @title Maximum eigen value normalisation
#'
#' @description Maximum eigen value row normalisation of a spatial weight matrix.
#
#' @param W spatial weights matrix
#'
#' @return
#' \describe{\emph{W}}  Eigen normalised spatial weights matrix
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{rownor}}
#'
#' @examples
#' data(gN3dist)
#' W<-InvDistMat(distMat=gN3dist,100000,powr=2,metres=FALSE)
#' Wnor<-eignor(W)
#'
#' @import Matrix
#' @import RSpectra
#'
#' @export


eignor<-function(W){
  if(nrow(W)!=ncol(W)) stop("Error in matrix!")
  Ws<-Matrix::Matrix(W, sparse = TRUE)
  emax<-Re(RSpectra::eigs(Ws,1,which = "LA")$values)
  W=W/emax
  return(W)
}
