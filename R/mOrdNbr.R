#' @name mOrdNbr
#' @title 1st to m-th order neighbours matrix
#'
#' @description Finds the from 1th to m-th order neighbours matrix.
#'
#' @param sf_pol spatial polygons object
#' @param m the order of neighbours up to which they will be included in the weights matrix, default 1
#' @param neigbs neighbours list, default NULL
#' @param listv logocal, default FALSE. If TRUE the list of neighbours should also be returned
#' @param rn logical, default FALSE. If TRUE, the weigth matrix will be row normalised
#' @param zrow logical, default TRUE. If FALSE, the row matrix is normalised even if there are zero rows in the matrix.
#'
#' @return
#' \describe{\emph{W}} spatial weights matrix (and list of neighbours \emph{nlist})
#'
#' @author Rozeta Simonovska
#'
#' @import spdep
#' @import methods
#'
#' @examples
#' library("rgdal")
#' ger<-readOGR(system.file(dsn="shape",package="SDPDmod"),layer="GermanyNUTS3")
#' m1thn<-mOrdNbr(ger)
#' m4thn<-mOrdNbr(ger,4)
#' mat1<-rownor(m4thn)
#' m4thn2<-mOrdNbr(ger,4,listv=TRUE,rn=TRUE)
#' mat2<-m4thn2$W
#'
#' @export

mOrdNbr<-function(sf_pol, m = 1, neigbs = NULL, listv = FALSE, rn = FALSE, zrow = TRUE){
  if(!is.list(neigbs) & length(neigbs)==0) {
    if(!is(sf_pol,"SpatialPolygons")) {stop("Wrong entry! Value must be a spatial polygons object.")}
    neigbs<-spdep::poly2nb(sf_pol)
  }
  N<-length(neigbs)
  W<-matrix(0,nrow=N,ncol=N)
  for(i in 1:N){  if(all(unlist(neigbs[[i]])!=0)){  W[i,unlist(neigbs[[i]])]<-1  }  }
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
        if(all(unlist(mneigb[[i]]))!=0){
          for(j in unlist(mneigb[[i]])){   v.p<-c(v.p,unlist(neigbs[[j]]))  }
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

  ###row-normalisation
  if(rn){ W<-rownor(W, zrow)  }
  if(listv){ return(list(W=W,nlist=nbrL))
  }else{ return(W)}
}

