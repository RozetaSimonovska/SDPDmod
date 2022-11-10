#' @name SharedBMat
#' @title Shared boundary matrix
#'
#' @description This function calculates the shared boundary matrix
#'
#' @param sf_pol spatial polygons, spatial lines object or spatial data frame
#' @param rn logical, default FALSE. If TRUE, the spatial weights matrix is row-normalized
#'
#' @return
#' \item{W}{spatial weights matrix (length of shared boundary between spatial units)}
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' library("sf")
#' \donttest{
#' ger   <- st_read(system.file(dsn = "shape/GermanyNUTS3.shp", package = "SDPDmod"))
#' bav <- ger[which(substr(ger$NUTS_CODE,1,3)=="DE2"),] ## Bavaria districts
#' W     <- SharedBMat(bav)}
#'
#' @import spdep
#' @import sp
#' @import sf
#' @import methods
#'
#' @export

SharedBMat<-function(sf_pol, rn = FALSE){
  if((is(sf_pol,"data.frame") &
      all(sf::st_geometry_type(sf_pol) %in% c("POLYGON","MULTIPOLYGON","LINESTRING","MULTILINESTRING","GEOMETRYCOLLECTION")))){
    sf_pol2 <- as(sf_pol, "Spatial")
    sf_obj <- sf_pol
  }else if(is(sf_pol,"SpatialPolygons") || is(sf_pol,"SpatialLines")) {

    sf_obj <- st_as_sf(sf_pol)

    if(is(sf_pol,"SpatialPolygons")){
      sf_pol2<- sf_pol
    }else if(is(sf_pol,"SpatialLines")){
      sf_pol2 <- sp::SpatialPolygons(
        lapply(1:length(sf_pol),
               function(x) Polygons(lapply(coordinates(sf_pol)[[x]],
                                           function(y) Polygon(y)), as.character(x))))
    }

  }else{
    stop("Wrong data type! Data must be a spatial polygons, spatial lines object or data frame containing geometry.")
  }

  N         <- nrow(sf_obj)
  neigbs    <- spdep::poly2nb(sf_pol2)
  line_sf_l <- line_sf_l2 <- vector(mode = "list", length = N)
  W         <- matrix(0,nrow = N,ncol=N)
  for(i in 1:N){  line_sf_l[[i]] <- line_sf_l2[[i]] <- vector(mode = "list", length = N)}  ##line_sf_l <- line_sf_l2 <- lapply(1:N,function(x) rep(NA,N))
  for(i in 1:N){
    if(all(neigbs[[i]])!=0){
      ngb_v<-neigbs[[i]][neigbs[[i]]>i]
      if(length(ngb_v)>0){
        for(j in ngb_v){
          line_sf_l[[i]][[j]]  <- sf::st_intersection(
                                        st_cast(sf_obj$geometry[i],to = "MULTILINESTRING")
                                       ,st_cast(sf_obj$geometry[j],to = "MULTILINESTRING"))
          temp<-NULL
          if(inherits(line_sf_l[[i]][[j]],"sfc_GEOMETRYCOLLECTION")){
            temp <- st_collection_extract(line_sf_l[[i]][[j]],"LINESTRING")
          }else if(inherits(line_sf_l[[i]][[j]],"sfc_MULTIPOLYGON")){
            temp <- st_cast(line_sf_l[[i]][[j]],to = "MULTILINESTRING")
          }else if(inherits(line_sf_l[[i]][[j]],c("sfc_POINT","sfc_POLYGON"))){
            temp <- st_cast(line_sf_l[[i]][[j]],to = "LINESTRING")
          }else if(inherits(line_sf_l[[i]][[j]],c("sfc_MULTILINESTRING","sfc_LINESTRING"))){
            temp<-line_sf_l[[i]][[j]]}
          if(!inherits(line_sf_l[[i]][[j]],c("sfc_POINT","sfc_MULTIPOINT"))){
            line_sf_l2[[i]][[j]] <- as_Spatial(temp)
            if(is(line_sf_l2[[i]][[j]],"SpatialLines")){
              W[i,j] <- sum(sp::SpatialLinesLengths(line_sf_l2[[i]][[j]]))
              W[j,i] <- W[i,j]
            }
          }else{
            W[i,j] <- 0.1
            W[j,i] <- W[i,j]
          }
        }
      }
    }
  }
  if(rn){ W<-rownor(W) }
  return(W)
}
