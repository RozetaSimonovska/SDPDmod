#' @name ExpDistMat
#' @title Exponential distance matrix
#'
#' @description This function calculates the (negative) exponential distance matrix,
#' with a given cut off distance and a positive exponent value.
#'
#' @param distMat distance matrix
#' @param distCutOff cut off distance in metres.
#' If no value is specified, then exponential decay is calulated for all values.
#' @param expn positive exponent, default 0.01
#' @param metres logical, default TRUE. If FALSE, the distance is transformed into kilometres
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
#' data(gN3dist)
#' W1<-ExpDistMat(distMat=gN3dist)
#' W2<-ExpDistMat(distMat=gN3dist,100000,expn=0.001,metres=FALSE)
#' Wnor<-eignor(W2)
#'
#' @export


ExpDistMat<-function(distMat,distCutOff = NULL,expn = 0.01,metres = TRUE){
  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)
    if(!metres){ distMat<-distMat/1000; distCutOff<-distCutOff/1000}
    for(i in 1:(n-1)){
      if(!is.null(distCutOff)){
        for(j in which(distMat[i,]<distCutOff)){
          if(j>i){
              W[i,j]<-exp(-distMat[i,j]*expn)
              W[j,i]<-exp(-distMat[i,j]*expn)
          }
        }
      } else{
        for(j in (i+1):n){
            W[i,j]<-exp(-distMat[i,j]*expn)
            W[j,i]<-exp(-distMat[i,j]*expn)
        }
      }
    }
  } else { stop("Error in distMat! Not a distance matrix.")}
  return(W)
}
