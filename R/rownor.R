#' @name rownor
#' @title Row normalisation
#'
#' @description Row normalisation of a spatial weights matrix.
#'
#' @param W spatial weights matrix
#'
#' @return
#' \describe{\emph{W}}  row normalised spatial weights matrix
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{eignor}}
#'
#' @examples
#' library("rgdal")
#' ger<-readOGR(system.file(dsn="shape",package="SDPDmod"),layer="GermanyNUTS3")
#' W<-mOrdNbr(ger,3)
#' Wnor<-rownor(W)
#'
#' @export

rownor<-function(W){
  if(nrow(W)!=ncol(W)) stop("Error in matrix!")
  N<-nrow(W)
  if(all(rowSums(W)!=0)) { W<-W/rowSums(W)
  } else {
    message("Zero row in matrix.")
  }

  return(W)
}
