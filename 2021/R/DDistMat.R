#' @name DDistMat
#' @title Double-Power Distance Weights Matrix
#'
#' @description This function calculates the double-power distance matrix,
#' for a given distance cut-off and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff distance cut-off. Default = half of the maximal distance from the distance matrix.
#' @param powr power (positive exponent), default 2
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalisation is performed.
#'
#' @return
#' \describe{\emph{W}}  weights matrix (not normalised)
#'
#' @details
#' W is an \emph{nxn} matrix with elements \eqn{w_{ij}}, \eqn{i, j = 1, ... n}, where
#' \deqn{w_{ij} = (1-(\frac{d_{ij}}{D})^p)^p, if 0 <= d_{ij} < D}
#' \deqn{w_{ij} = 0, if d_{ij} > D or i = j}
#' where D is the cut-off distance point (maximum radius of influence),
#' \eqn{d_{ij}} is the distance between spatial units \emph{i} and \emph{j},
#' and \emph{p} is the power value (e.g. \emph{p} = 2, 3, 4,...).
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' data(gN3dist) ##distance in metres
#' W1<-DDistMat(distMat=gN3dist)
#' dist2<-gN3dist/1000 ##in km
#' W2<-DDistMat(distMat=gN3dist,300000,powr=3) ##distance in metres
#' W3<-DDistMat(distMat=dist2,300,powr=3)  ##distance in kilometres
#'
#' @export



DDistMat<-function(distMat, distCutOff = NULL, powr = 2, mevn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    if(is.null(distCutOff)){distCutOff <- max(distMat)/2 }
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)

    for(i in 1:(n-1)){
      for(j in which(distMat[i,]<distCutOff)){
        if(j>i){
          temp <- (1-(distMat[i,j]/distCutOff)^powr)^powr
          W[i,j]<-temp
          W[j,i]<-temp
        }
      }
    }

    if(mevn){ W <- eignor(W) }

  } else { stop("Error in distMat! Not a distance matrix.")}

  return(W)
}
