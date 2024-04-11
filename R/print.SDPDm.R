#' @name print.SDPDm
#'
#' @title print for class SDPDm
#'
#' @description Method for sprinting the results of objects of class "SDPDm"
#'
#' @method print SDPDm
#'
#' @param x object of class "SDPDm"
#' @param digits number of digits
#' @param ... additional arguments to be passed
#'
#' @return No return value
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#' 
#' @export

print.SDPDm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if(inherits(x,"SDPDm")){
    if(x$dynamic){
      cat(paste0(x$model," dynamic panel model with ",
                 x$effect, " fixed effects\n"))
    }else{
      cat(paste0(x$model," panel model with ",
                 x$effect, " fixed effects\n"))
    }
    cat("\nCall:\n")
    print(x$call)
    cat("\nSpatial autoregressive coefficient:\n")
    
    if((x$dynamic & x$LeeYu & x$effect %in%
        c("individual","twoways")) ||
       (x$dynamic & x$DirectT & x$effect %in% c("twoways"))){
      
      cat(round(x$rho1,digits))
    }else{
      cat(round(x$rho,digits))
    }

    cat("\nCoefficients:\n")
    if((x$dynamic & x$LeeYu & x$effect %in%
        c("individual","twoways")) ||
       (x$dynamic & x$DirectT & x$effect %in% c("twoways"))){
      print(round(x$coefficients1,digits))
    }else{
      print(round(x$coefficients,digits))
    }
    cat("\n")
  }
  invisible(x)
}
