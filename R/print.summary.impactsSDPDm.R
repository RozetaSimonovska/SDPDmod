#' @name print.summary.impactsSDPDm
#'
#' @title Print summary for class impactsSDPDm
#'
#' @description Method for printing the summary the results of objects of class "impactsSDPDm"
#'
#' @method print summary.impactsSDPDm
#'
#' @param x summary object of class "impactsSDPDm"
#' @param ... additional arguments to be passed
#'
#' @author Rozeta Simonovska
#'
#' @importFrom stats symnum
#' 
#' @export

print.summary.impactsSDPDm <- function(x,...) {

    if(length(x)==14){
      cat("\nImpact estimates for spatial dynamic model\n")
      cat("========================================================\n")
      cat("Short-term\n")
      cat("\nDirect:\n")
      printCoefmat(x$DIRECTst.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(x$INDIRECTst.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(x$TOTALst.tab, signif.legend = FALSE)
      cat("========================================================\n")
      cat("Long-term\n")
      cat("\nDirect:\n")
      printCoefmat(x$DIRECTlt.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(x$INDIRECTlt.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(x$TOTALlt.tab, signif.legend = FALSE)
      
      p.values <- c(x$DIRECTst.tab[,4],x$INDIRECTst.tab[,4],
                    x$TOTALst.tab[,4],
                    x$DIRECTlt.tab[,4],x$INDIRECTlt.tab[,4],
                    x$TOTALlt.tab[,4]
      )
      if(length(p.values[which(p.values<0.1)])>0){
        Signif <- stats::symnum(p.values, corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
        sl <- attr(Signif, "legend")
        cat("---\nSignif. codes:  ", sl, sep = "",
            fill = getOption("width") + 4 + max(nchar(sl, "bytes") - nchar(sl)))
      }
    }else if(length(x)==7){
      cat("\nImpact estimates for spatial (static) model\n")
      cat("========================================================\n")
      cat("\nDirect:\n")
      printCoefmat(x$DIRECT.tab, signif.legend = FALSE)
      cat("\nIndirect:\n")
      printCoefmat(x$INDIRECT.tab, signif.legend = FALSE)
      cat("\nTotal:\n")
      printCoefmat(x$TOTAL.tab, signif.legend = FALSE)
      
      p.values <- c(x$DIRECT.tab[,4],x$INDIRECT.tab[,4],
                    x$TOTAL.tab[,4]
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
  
  invisible(x)
}
