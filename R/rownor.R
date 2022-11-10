#' @name rownor
#' @title Row-normalization
#'
#' @description Row-normalization of a spatial weights matrix.
#'
#' @param W spatial weights matrix
#'
#' @return
#' \item{W}{row-normalized spatial weights matrix}
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{eignor}}
#'
#' @examples
#' library("sf")
#' ger   <- st_read(system.file(dsn = "shape/GermanyNUTS3.shp",package = "SDPDmod"))
#' W     <- mOrdNbr(ger, 3)
#' Wnor  <- rownor(W)
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
