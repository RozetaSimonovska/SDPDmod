#' @name coef.SDPDm
#'
#' @title Extract coefficents from model of class SDPDm
#'
#' @description Method for extracting coefficients of objects of class "SDPDm"
#'
#' @method coef SDPDm
#'
#' @param object object of class "SDPDm"
#' @param ... additional arguments to be passed
#'
#' @return Coefficients extracted from the model object of class "SDPDm".
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#'
#' @export


coef.SDPDm <- function(object,...) {
  if(inherits(object,"SDPDm")){
    if((object$dynamic & object$LeeYu & object$effect %in%
        c("individual","twoways")) ||
       (object$dynamic & object$DirectT & object$effect %in% c("twoways"))){
      return(object$coefficients1)
    }else{
      return(object$coefficients)
    }
  }
}
