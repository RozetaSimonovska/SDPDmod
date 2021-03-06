#' @name ExpDistMat
#' @title Exponential distance matrix
#'
#' @description This function calculates the (negative) exponential distance matrix,
#' with a given cut-off distance and a positive exponent value.
#'
#' @param distMat distance matrix
#' @param distCutOff cut-off distance. Default = half of the maximal distance from the distance matrix.
#' @param expn positive exponent, default = 0.01
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalisation is performed.
#'
#' @return
#' \describe{\emph{W}}  weights matrix (not normalised)
#'
#' @details
#' W is an \emph{nxn} matrix with elements \eqn{w_{ij}}, \emph{i, j = 1,..n}, where
#' \deqn{w_{ij}=e^{-\alpha d_{ij}}, if 0 <= d_{ij} < D}
#' \deqn{w_{ij}=0, if d_{ij} > D or i = j}
#' where D is the distance cut-off point (maximum radius of influence),
#' \eqn{d_{ij}} is the distance between spatial units \emph{i} and \emph{j}, and
#' \eqn{\alpha} is the positive exponent (e.g. \eqn{\alpha}= 0.01, 0.02,...).
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' data(gN3dist) ##distance in metres
#' W1<-ExpDistMat(distMat=gN3dist, distCutOff=100000)
#' dist2<-gN3dist/1000 ##in km
#' W2<-ExpDistMat(distMat=dist2, distCutOff=100, expn=0.001)
#' W2nor<-ExpDistMat(distMat=dist2, distCutOff=100000, expn=0.001, mevn=TRUE)
#'
#' @export


ExpDistMat<-function(distMat, distCutOff = NULL, expn = 0.01, mevn = FALSE){
  if(isSymmetric(distMat) & all(diag(distMat)==0)){

    if(is.null(distCutOff)){distCutOff <- max(distMat)/2 }
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)

    for(i in 1:(n-1)){
      for(j in which(distMat[i,]<distCutOff & distMat[i,]!=0)){
          if(j>i){
              temp <- exp(-distMat[i,j]*expn)
              W[i,j]<- temp
              W[j,i]<- temp
          }
      }
    }

    if(mevn){ W <- eignor(W) }

  } else { stop("Error in distMat! Not a distance matrix.")}

  return(W)
}
