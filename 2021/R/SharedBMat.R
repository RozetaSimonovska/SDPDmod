#' @name SharedBMat
#' @title Shared boundary matrix
#'
#' @description This function calculates the shared boundary matrix
#'
#' @param sf_pol spatial polygons or spatial lines object
#' @param rn logical, default FALSE. If TRUE, the weigth matrix is row normalised
#'
#' @return
#' \describe{\emph{W}} spatial weights matrix (length of shared boundary between spatial units)
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' library("rgdal")
#' ger<-readOGR(system.file(dsn="shape",package="SDPDmod"),layer="GermanyNUTS3")
#' W<-SharedBMat(ger)
#'
#' @import spdep
#' @import sp
#' @import rgeos
#' @import methods
#'
#' @export

SharedBMat<-function(sf_pol, rn = FALSE){

  if(!is(sf_pol,"SpatialPolygons") & !is(sf_pol,"SpatialLines")) {
    stop("Wrong entry! Value must be a spatial polygons or spatial lines object.")
  }else if(is(sf_pol,"SpatialPolygons")) {
    sf_line <- as(sf_pol, "SpatialLines")
    sf_pol2<-sf_pol
  }else {
    sf_line <- sf_pol
    sf_pol2 <- sp::SpatialPolygons(
      lapply(1:length(sf_line),
             function(x) Polygons(lapply(coordinates(sf_line)[[x]], function(y) Polygon(y)), as.character(x))))
    }

  N<-length(sf_line)
  neigbs<-spdep::poly2nb(sf_pol2)
  line_sf_l<-vector(mode = "list", length = N)
  W<-matrix(0,nrow = N,ncol=N)
  for(i in 1:N){line_sf_l[[i]]<-vector(mode = "list", length = N)}
  for( i in 1:N){
    if(all(neigbs[[i]])!=0){
      for(j in neigbs[[i]]){
        line_sf_l[[i]][[j]]<- rgeos::gIntersection(sf_line[i,], sf_line[j,],byid = TRUE)
        if(is(line_sf_l[[i]][[j]],"SpatialLines")){
          W[i,j] <- sp::SpatialLinesLengths(line_sf_l[[i]][[j]])  }
      }
    }
  }
  if(rn){ W<-rownor(W) }
  return(W)
}
