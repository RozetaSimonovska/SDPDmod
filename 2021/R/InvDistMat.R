#' @name InvDistMat
#' @title Inverse distance matrix
#'
#' @description This function calculates the inverse distances,
#' with a given cut-off distance and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff cut-off distance. Default = half of the maximal distance from the distance matrix.
#' @param powr power (positive exponent), default = 1
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalisation is performed.
#'
#' @return
#' \describe{\emph{W}}  weights matrix (not normalised)
#'
#' @details
#' W is an \emph{nxn} matrix with elements \eqn{w_{ij}}, \emph{i,j=1,..n}, where
#' \deqn{w_{ij}=1/d_{ij}^\gamma, if 0 <= d_{ij} < D}
#' \deqn{w_{ij}=0, if d_{ij} > D or i = j}
#' where \emph{D} is the distance cut-off point (maximum radius of influence),
#' \eqn{d_{ij}} is the distance between spatial units \emph{i} and \emph{j},
#' and \eqn{\gamma} is the value for the exponent (e.g. \eqn{\gamma} = 1, 2, 3, 4,...).
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' ## distance between centroids of NUTS3 regions in Germany (in metres)
#' data(gN3dist,package = "SDPDmod")
#' ## inverse distance matrix with cut-off 100000 metres
#' W1<-InvDistMat(distMat=gN3dist, distCutOff=100000)
#' dist2<-gN3dist/1000 ##distance in km
#' ## normalised distance matrix with cut-off 100km
#' W2<-InvDistMat(distMat=dist2, distCutOff=100, powr=2, mevn=TRUE)
#'
#' @export


InvDistMat<-function(distMat, distCutOff = NULL, powr = 1, mevn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    if(is.null(distCutOff)){distCutOff <- max(distMat)/2 }
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)

    for(i in 1:(n-1)){
      for(j in which(distMat[i,]<distCutOff & distMat[i,]!=0)){
          if(j>i){
            temp <- 1/(distMat[i,j])^powr
            W[i,j]<- temp
            W[j,i]<- temp
          }
      }
    }

    if(mevn){ W <- eignor(W) }

  } else { stop("Error in distMat! Not a distance matrix.")}

  return(W)
}
