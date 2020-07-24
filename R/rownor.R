#' @name rownor
#' @title Row normalisation
#'
#' @description Row normalisation of a spatial weights matrix.
#'
#' @param W spatial weights matrix
#' @param zrow logical, default TRUE. If FALSE, the row matrix is normalised even if there are zero rows in the matrix.
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

rownor<-function(W,zrow = TRUE){
  if(nrow(W)!=ncol(W)) stop("Error in matrix!")
  N<-nrow(W)
  if(all(rowSums(W)!=0)) { W<-W/rowSums(W)
  } else {
    if(zrow){message("Zero row in matrix.")
    }else{
      for(i in 1:N){ if(sum(W[i,])!=0) W[i,]=W[i,]/sum(W[i,])  }
    }
  }
  return(W)
}
