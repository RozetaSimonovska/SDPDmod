#' @name print.blmpSDPD
#'
#' @title Print for class blmpSDPD
#'
#' @description Method for printing the results of objects of class "blmpSDPD"
#'
#' @method print blmpSDPD
#'
#' @param x object of class "blmpSDPD"
#' @param digits number of digits
#' @param ... additional arguments to be passed
#'
#' @return No return value
#' 
#' @author Rozeta Simonovska
#' 
#' @export

print.blmpSDPD <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    
    if(inherits(x,"blmpSDPD")){
      res <- rbind(x$lmarginal,x$probs)
      rownames(res) <- c("Log marginals","Model probs")
      
      print(round(res,digits))
    }
    invisible(x)
  }