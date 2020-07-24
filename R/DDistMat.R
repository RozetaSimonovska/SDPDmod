#' @name DDistMat
#' @title Double-Power Distance Weights Matrix
#'
#' @description This function calculates the double-power distance matrix,
#' for a given cut off distance and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff cut off distance (in metres). Default 100000 metres
#' @param powr power (positive exponent), default 2
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
#' data(gN3dist)
#' W1<-DDistMat(distMat=gN3dist)
#' dist2<-gN3dist/1000 ##in km
#' W2<-DDistMat(distMat=gN3dist,300000,powr=3) ##distance in metres
#' W3<-DDistMat(distMat=dist2,300,powr=3)  ##distance in kilometres
#'
#' @export



DDistMat<-function(distMat,distCutOff = 100000,powr = 2){
  if(isSymmetric(distMat) & all(diag(distMat)==0)){
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
  } else { stop("Error in distMat! Not a distance matrix.")}
  return(W)
}
