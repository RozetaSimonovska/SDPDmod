#' @name DDistMat
#' @title Double-Power Distance Weights Matrix
#'
#' @description This function calculates the double-power distance matrix,
#' for a given distance cutoff and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff distance cutoff. Default = the maximal value from the distance matrix.
#' @param powr power (positive exponent), default 2
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalization is performed.
#'
#' @return
#' \item{W}{spatial weights matrix (Default, not normalized)}
#'
#' @details
#' W is an \emph{nxn} matrix with elements \eqn{w_{ij}}, \eqn{i, j = 1, ... n}, where
#' \eqn{w_{ij} = (1-(\frac{d_{ij}}{D})^p)^p}, if \eqn{0 <= d_{ij} < D} and
#' \eqn{w_{ij} = 0}, if \eqn{d_{ij} > D} or \eqn{i = j}.
#' \emph{D} is the cut-off distance point (maximum radius of influence),
#' \eqn{d_{ij}} is the distance between spatial units \emph{i} and \emph{j},
#' and \emph{p} is the power value (e.g. \emph{p} = 2, 3, 4,...).
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' data(gN3dist) ##distance in meters
#' W1    <- DDistMat(distMat = gN3dist, distCutOff = 300000, powr = 3) ##distance cutoff in meters
#' dist2 <- gN3dist/1000 ##in km
#' W2    <- DDistMat(distMat = dist2, 300, 3)  ##distance cutoff in kilometers
#'
#' @export



DDistMat<-function(distMat, distCutOff = NULL, powr = 2, mevn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    if(is.null(distCutOff)){distCutOff <- max(distMat) }
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
