#' @name summary.impactsSDPDm
#'
#' @title Summary for class impactsSDPDm
#'
#' @description Method for summarizing the results of objects of class "impactsSDPDm"
#'
#' @method summary impactsSDPDm
#'
#' @param object object of class "impactsSDPDm"
#' @param ... additional arguments to be passed
#'
#' @return No return value
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#'
#' @importFrom stats symnum
#'
#' @export


summary.impactsSDPDm <- function(object,...) {
  if(inherits(object,"impactsSDPDm")){
    if(length(object)==6){
      cat("\nImpact estimates for spatial dynamic model\n")
      cat("========================================================\n")
      cat("Short-term\n")
      cat("\nDirect:\n")
      printCoefmat(object$DIRECTst.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(object$INDIRECTst.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(object$TOTALst.tab, signif.legend = FALSE)
      cat("========================================================\n")
      cat("Long-term\n")
      cat("\nDirect:\n")
      printCoefmat(object$DIRECTlt.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(object$INDIRECTlt.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(object$TOTALlt.tab, signif.legend = FALSE)

      p.values <- c(object$DIRECTst.tab[,4],object$INDIRECTst.tab[,4],
                    object$TOTALst.tab[,4],
                    object$DIRECTlt.tab[,4],object$INDIRECTlt.tab[,4],
                    object$TOTALlt.tab[,4]
                    )
      if(length(p.values[which(p.values<0.1)])>0){
        Signif <- stats::symnum(p.values, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
        sl <- attr(Signif, "legend")
        cat("---\nSignif. codes:  ", sl, sep = "",
            fill = getOption("width") + 4 + max(nchar(sl, "bytes") - nchar(sl)))
      }
    }else if(length(object)==3){
      cat("\nImpact estimates for spatial (static) model\n")
      cat("========================================================\n")
      cat("\nDirect:\n")
      printCoefmat(object$DIRECT.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(object$INDIRECT.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(object$TOTAL.tab, signif.legend = FALSE)

      p.values <- c(object$DIRECT.tab[,4],object$INDIRECT.tab[,4],
                    object$TOTAL.tab[,4]
                    )
      if(length(p.values[which(p.values<0.1)])>0){
        Signif <- stats::symnum(p.values, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
        sl <- attr(Signif, "legend")
        cat("---\nSignif. codes:  ", sl, sep = "",
            fill = getOption("width") + 4 + max(nchar(sl, "bytes") - nchar(sl)))
      }

    }
  }
}
