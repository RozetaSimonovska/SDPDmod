#' @name mOrdNbr
#' @title 1st to m-th order neighbors matrix
#'
#' @description Finds the 1th to m-th order neighbors matrix.
#'
#' @param sf_pol spatial polygons object
#' @param m the order of neighbors up to which they will be included in the weights matrix, default 1
#' @param neigbs neighbors list, default NULL
#' @param listv logical, default FALSE. If TRUE the list of neighbors should also be returned
#' @param rn logical, default FALSE. If TRUE, the weight matrix will be row-normalized
#'
#' @return
#' \item{W}{spatial weights matrix}
#' \item{nlist}{list of neighbors}
#'
#' @author Rozeta Simonovska
#'
#' @import sf
#' @import spdep
#' @import methods
#'
#' @examples
#' library("sf")
#' ger   <- st_read(system.file(dsn = "shape/GermanyNUTS3.shp",package = "SDPDmod"))
#' m1thn <- mOrdNbr(ger)
#' \donttest{m4thn <- mOrdNbr(ger, 4)
#' mat1  <- rownor(m4thn)
#' m4thn2<- mOrdNbr(ger, 4, listv = TRUE, rn = TRUE)
#' mat2  <- m4thn2$W}
#'
#' @export

mOrdNbr<-function(sf_pol=NULL, m = 1, neigbs = NULL, listv = FALSE, rn = FALSE){

  if(is.null(sf_pol) & is.null(neigbs)){
    stop("Missing value for sf_pol and neigbs! At least one value for sf_pol or neigbs has to be entered.")
  }else if(!is.null(neigbs) & !is.list(neigbs) & length(neigbs)==0) {
    stop("Error in neighbours")
  }else if(!is.null(sf_pol)){
    if(!is(sf_pol,"SpatialPolygons") &
         !(is(sf_pol,"data.frame") & all(sf::st_geometry_type(sf_pol) %in% c("POLYGON","MULTIPOLYGON")))) {
        stop("Wrong data type! Data must be a spatial polygons object or data frame containing geometry.")
    }else{

    if(is(sf_pol,"data.frame") & all(sf::st_geometry_type(sf_pol) %in% c("POLYGON","MULTIPOLYGON"))){
        sf_pol2 <- as(sf_pol, "Spatial")
      }else{ sf_pol2 <- sf_pol}

      neigbs<-spdep::poly2nb(sf_pol2)
    }
  }


  N<-length(neigbs)
  W<-matrix(0,nrow=N,ncol=N)

  for(i in 1:N){  if(all(neigbs[[i]]!=0)){  W[i,neigbs[[i]]]<-1  }  }
  nbrL<-vector("list",m)
  nbrL[[1]]<-neigbs
  if(m>1){
    for(j in 2:m){ nbrL[[j]]<-vector("list",N) }
    k<-2
    repeat{
      for(i in 1:N){
        v.p<-vector()
        mneigb<-nbrL[[k-1]]
        v.n<-as.list(1:N)
        if(all(mneigb[[i]])!=0){
          for(j in mneigb[[i]]){   v.p<-c(v.p,neigbs[[j]])  }
          v.pp<-unique(v.p)
          v.pp<-v.pp[order(v.pp)]
          for(l in 1:(m-1)){  v.n[[i]]<-c(v.n[[i]],nbrL[[l]][[i]])    }
          v.ppp<-v.pp[which(!v.pp %in% c(i,v.n[[i]]))]
          nbrL[[k]][[i]]<-v.ppp
          if(length(v.ppp)!=0){   W[i,v.ppp]<-1   }
        }
      }
      k<-k+1
      if(k>m){break}
    }
  }

  ###row-normalization
  if(rn){ W<-rownor(W)  }

  if(listv){
    return(list(W=W,nlist=nbrL))
  }else{
    return(W)
    }
}

