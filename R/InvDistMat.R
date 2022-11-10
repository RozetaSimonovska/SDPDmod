#' @name InvDistMat
#' @title Inverse distance matrix
#'
#' @description This function calculates the inverse distances,
#' with a given cutoff distance and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff cutoff distance. Default = the maximal value from the distance matrix.
#' @param powr power (positive exponent), default = 1
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalization is performed.
#'
#' @return
#' \item{W}{weights matrix (Default, not normalized)}
#'
#' @details
#' W is an \emph{nxn} matrix with elements \eqn{w_{ij}}, \emph{i,j=1,..n}, where
#' \eqn{w_{ij}=1/d_{ij}^\gamma}, if \eqn{0 <= d_{ij} < D} and
#' \eqn{w_{ij}=0}, if \eqn{d_{ij} > D} or \eqn{i = j}.
#' \emph{D} is the distance cutoff point (maximum radius of influence),
#' \eqn{d_{ij}} is the distance between spatial units \emph{i} and \emph{j},
#' and \eqn{\gamma} is the value for the exponent (e.g. \eqn{\gamma} = 1, 2, 3, 4,...).
#'
#' @author Rozeta Simonovska
#'
#' @examples
#' ## distance between centroids of NUTS3 regions in Germany (in meters)
#' data(gN3dist,package = "SDPDmod")
#' ## inverse distance matrix with cutoff 100000 meters
#' W1    <- InvDistMat(distMat = gN3dist, distCutOff = 100000)
#' dist2 <- gN3dist/1000 ##distance in km
#' ## normalized distance matrix with cutoff 100km
#' W2    <- InvDistMat(distMat = dist2, distCutOff=100, powr = 2, mevn = TRUE)
#'
#' @export


InvDistMat<-function(distMat, distCutOff = NULL, powr = 1, mevn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    if(is.null(distCutOff)){distCutOff <- max(distMat) }
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
