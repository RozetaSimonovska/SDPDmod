#' @name mNearestN
#'
#' @title m nearest neigbours based on a distance matrix
#'
#' @description This function finds the m nearest neaigbours, given a matrix of distances.
#'
#' @param distMat distance matrix
#' @param m number of nearest neigbours, default value 5
#' @param listv logocal, default FALSE. If TRUE the list of neighbours should also be returned
#' @param rn logical, default FALSE. If TRUE, the weigth matrix will be row normalised
#'
#' @return
#' \describe{\emph{W}} weights matrix
#' \describe{\emph{nlist}} list of indexes of the m nearest neigbours
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' data(gN3dist,package = "SDPDmod")
#' fournn<-mNearestN(gN3dist,4)
#' mat1<-rownor(fournn)
#' tennn<-mNearestN(gN3dist,10,listv=TRUE,rn=TRUE)
#' mat2<-tennn$W
#'
#' @export

###Funcrion

mNearestN<-function(distMat, m = 5,listv = FALSE, rn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    n<-nrow(distMat)
    list_ngb<-vector("list",n)

    for(i in 1:n){
      ordRow<-order(distMat[i,], decreasing = FALSE)
      nRow<-ordRow[which(ordRow!=i)]
      if(length(nRow)>=m){
        list_ngb[[i]]<-nRow[1:m]
      } else {
        list_ngb[[i]]<-nRow[1:n]
        warning("m is larger than total number of individuals/regions")
      }
    }
  } else { stop("Error in distMat! Not a distance matrix.")}

  W<-matrix(0,nrow=n,ncol=n)
  for(i in 1:n){  W[i,list_ngb[[i]]]<-1  }
  if(rn){ W <- rownor(W) }

  if(listv){ return(list(W = W, nlist = list_ngb))
    }else{ return(W) }
}
####END of function
