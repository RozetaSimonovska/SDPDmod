##beta prior centered on zero ##a0=c(1.01,1.1,2)
betapr<-function(x, a0 = 1.01){
  if(!is.numeric(a0)){ message("Error in a0! Non-numeric value.")
  }else{
    prior<-(1/beta(a0,a0))*((((1-x)^(a0-1))*((1+x)^(a0-1)))/(2^(a0+a0-1)))
    return(prior)
  }
}
