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
#' dist2<-gN3dist/1000 ##distance in km
#' W<-InvDistMat(distMat = dist2, distCutOff = 100, powr = 2)
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
