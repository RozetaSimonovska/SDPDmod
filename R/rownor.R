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
#' ger   <- st_read(system.file("shape/GermanyNUTS3.shp",
#'                              package = "SDPDmod"),
#'                  quiet = TRUE)
#' W     <- mOrdNbr(ger, 3)
#' Wnor  <- rownor(W)
#'
#' @export

rownor<-function(W){
  if(nrow(W)!=ncol(W)) stop("Error in matrix!")
  N <- nrow(W)
  Wn <- matrix(NA,N,N)
  if(all(rowSums(W)!=0)) {
    Wn<-W/rowSums(W)
  } else {
    warning("Zero row in matrix.")
    for(i in 1:N){
      if(sum(W[i,])!=0) {
        Wn[i,]<-W[i,]/sum(W[i,])
      }else{
        Wn[i,]<-0
      }
    }
  }

  return(Wn)
}
