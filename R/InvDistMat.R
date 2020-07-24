#' @name InvDistMat
#' @title Inverse distance matrix
#'
#' @description This function calculates the inverse distances,
#' with a given cut off distance and a positive exponent.
#'
#' @param distMat distance matrix
#' @param distCutOff cut off distance (in metres)
#' @param powr power (positive exponent), default 1
#' @param metres logical, default TRUE. If FALSE, the distance is transformed into kilometres
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
#' data(gN3dist,package = "SDPDmod")  ###distance between centroids of NUTS3 regions in Germany
#' W1<-InvDistMat(distMat=gN3dist)
#' W2<-InvDistMat(distMat=gN3dist,100000,powr=2,metres=FALSE)
#'
#' @export


InvDistMat<-function(distMat,distCutOff = NULL,powr = 1,metres = TRUE){
  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)
    if(!metres){ distMat<-distMat/1000; distCutOff<-distCutOff/1000}
    for(i in 1:(n-1)){
      if(!is.null(distCutOff)){
        for(j in which(distMat[i,]<distCutOff)){
          if(j>i){
            W[i,j]<-1/(distMat[i,j])^powr
            W[j,i]<-1/(distMat[i,j])^powr
          }
        }
      } else{
        for(j in (i+1):n){
          W[i,j]<-1/(distMat[i,j])^powr
          W[j,i]<-1/(distMat[i,j])^powr
        }
      }
    }
  } else { stop("Error in distMat! Not a distance matrix.")}
  return(W)
}
