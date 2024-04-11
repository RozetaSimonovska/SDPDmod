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
#' @return Summary of impacts
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#'
#' @export


summary.impactsSDPDm <- function(object,...) {
  if(inherits(object,"impactsSDPDm")){
    if(length(object)==6){
     
      est.st <- as.data.frame(matrix(c(object$DIRECTst.tab[,1],
                                       object$INDIRECTst.tab[,1],
                                       object$TOTALst.tab[,1]), ncol = 3))
      est.lt <- as.data.frame(matrix(c(object$DIRECTlt.tab[,1],
                                       object$INDIRECTlt.tab[,1],
                                       object$TOTALlt.tab[,1]), ncol = 3))
      
      sterr.st <- as.data.frame(matrix(c(object$DIRECTst.tab[,2],
                                       object$INDIRECTst.tab[,2],
                                       object$TOTALst.tab[,2]), ncol = 3))
      sterr.lt <- as.data.frame(matrix(c(object$DIRECTlt.tab[,2],
                                       object$INDIRECTlt.tab[,2],
                                       object$TOTALlt.tab[,2]), ncol = 3))
      
      tst.st <- as.data.frame(matrix(c(object$DIRECTst.tab[,3],
                                       object$INDIRECTst.tab[,3],
                                       object$TOTALst.tab[,3]), ncol = 3))
      tst.lt <- as.data.frame(matrix(c(object$DIRECTlt.tab[,3],
                                       object$INDIRECTlt.tab[,3],
                                       object$TOTALlt.tab[,3]), ncol = 3))
      
      p.values.st <- as.data.frame(matrix(c(object$DIRECTst.tab[,4],
                                       object$INDIRECTst.tab[,4],
                                       object$TOTALst.tab[,4]), ncol = 3))
      p.values.lt <- as.data.frame(matrix(c(object$DIRECTlt.tab[,4],
                                       object$INDIRECTlt.tab[,4],
                                       object$TOTALlt.tab[,4]), ncol = 3))
      
      rownames(est.st) <- rownames(sterr.st) <- 
        rownames(tst.st) <- rownames(p.values.st)<-
        rownames(est.lt) <- rownames(sterr.lt) <- 
        rownames(tst.lt) <- rownames(p.values.lt)<-
        rownames(object$TOTALst.tab)
      colnames(est.st) <- colnames(sterr.st) <- 
        colnames(tst.st) <- colnames(p.values.st)<-
        colnames(est.lt) <- colnames(sterr.lt) <- 
        colnames(tst.lt) <- colnames(p.values.lt)<-
        c("Direct","Indirect","Total")
      
      object$estimates.st <- est.st
      object$sterrors.st <-sterr.st
      object$tstats.st <- tst.st
      object$pvalues.st <- p.values.st
      object$estimates.lt <- est.lt
      object$sterrors.lt <-sterr.lt
      object$tstats.lt <- tst.lt
      object$pvalues.lt <- p.values.lt
  
    }else if(length(object)==3){
      
      est <- as.data.frame(matrix(c(object$DIRECT.tab[,1],object$INDIRECT.tab[,1],
                      object$TOTAL.tab[,1]), ncol = 3))
      
      sterr <- as.data.frame(matrix(c(object$DIRECT.tab[,2],object$INDIRECT.tab[,2],
                      object$TOTAL.tab[,2]), ncol = 3))
      
      tst <- as.data.frame(matrix(c(object$DIRECT.tab[,3],object$INDIRECT.tab[,3],
                      object$TOTAL.tab[,3]), ncol = 3))

      p.values <- as.data.frame(matrix(c(object$DIRECT.tab[,4],object$INDIRECT.tab[,4],
                    object$TOTAL.tab[,4]), ncol = 3))
      
      rownames(est) <- rownames(sterr) <- rownames(tst) <- rownames(p.values)<-
        rownames(object$TOTAL.tab)
      colnames(est) <- colnames(sterr) <- colnames(tst) <- colnames(p.values)<-
        c("Direct","Indirect","Total")
      
      object$estimates <- est
      object$sterrors <-sterr
      object$tstats <- tst
      object$pvalues <- p.values
    }
    class(object) <- c("summary.impactsSDPDm","impactsSDPDm")
    object
  }

}
