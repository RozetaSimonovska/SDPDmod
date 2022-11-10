#' @name DistWMat
#' @title Distance weights matrix (Inverse distance, Exponential distance or Double-Distance matrix)
#'
#' @description This function calculates the spatial distance weights matrix (inverse, exponential or
#' double-distance), with a given cutoff distance and a positive exponent (alpha).
#'
#' @param distMat distance matrix
#' @param distCutOff cutoff distance. Default = the maximal value from the distance matrix.
#' @param type the type of distance matrix c("inverse","expo","doubled"). Default = "inverse".
#' @param alpha power (positive exponent), default 1 if type="inverse", 0.01 if type="expo" and 2 if type="double"
#' @param mevn logical, default FALSE. If TRUE, max-eigenvalue normalization is performed.
#'
#' @return
#' \item{W}{spatial weights matrix (Default, not normalized)}
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{InvDistMat}} \code{\link{ExpDistMat}} \code{\link{DDistMat}} \code{vignette("spatial_matrices", package = "SDPDmod")}
#'
#' @examples
#' ## distance between centroids of NUTS3 regions in Germany (in meters)
#' data(gN3dist,package = "SDPDmod")
#' ##inverse distance matrix with cutoff 100000 meters
#' W1    <- DistWMat(distMat = gN3dist, distCutOff = 100000)
#' dist2 <- gN3dist/1000 ##distance in km
#' ## normalized exponential distance matrix
#' W2    <- DistWMat(distMat=dist2, distCutOff = 100, type = "expo", alpha = 2, mevn = TRUE)
#'
#' @export


DistWMat<-function(distMat, distCutOff = NULL, type= "inverse", alpha = NULL, mevn = FALSE){

  if(isSymmetric(distMat) & all(diag(distMat)==0)){
    if(is.null(distCutOff)){distCutOff <- max(distMat) }
    n<-nrow(distMat)
    W<-matrix(0,nrow = n,ncol = n)

    for(i in 1:(n-1)){
      for(j in which(distMat[i,]<distCutOff & distMat[i,]!=0)){
        if(j>i){
          if(type=="inverse"){
              if(is.null(alpha)){ alpha <- 1 }
              temp <- 1/(distMat[i,j])^alpha
          }else if(type=="expo"){
              if(is.null(alpha)){ alpha <- 0.01 }
              temp <- exp(-distMat[i,j]*alpha)
          }else if(type=="doubled"){
            if(is.null(alpha)){ alpha <- 2 }
            temp <- (1-(distMat[i,j]/distCutOff)^alpha)^alpha
          }else{  stop("Error in type of distance matrix!")  }
          W[i,j]<- temp
          W[j,i]<- temp
        }
      }
    }

    if(mevn){ W <- eignor(W) }

  } else { stop("Error in distMat! Not a distance matrix.")}

  return(W)
}
