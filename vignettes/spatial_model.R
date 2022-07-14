## ----setup, include = FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----------------------------------------------------------------------------------------------------
library(SDPDmod)

## ----data, eval=TRUE, echo=T, warning=FALSE, message=FALSE-------------------------------------------
data("Cigar",package = "plm")
head(Cigar)
data1<- Cigar
data1$logc<-log(data1$sales)
data1$logp<-log(data1$price/data1$cpi)
data1$logy<-log(data1$ndi/data1$cpi)
data1$lpm<-log(data1$pimin/data1$cpi)

data("usa46",package="SDPDmod") ## binary contiguity matrix of 46 USA states
str(usa46)
W <- rownor(usa46) ## row-normalization
isrownor(W) ## check if W is row-normalized

## ----m1, eval=TRUE, echo=T, warning=FALSE, message=FALSE---------------------------------------------
res1<-blmpSDPD(formula = logc ~ logp+logy, data = data1, W = W,
               index = c("state","year"),
               model = list("ols","sar","sdm","sem","sdem","slx"), 
               effect = "individual")
res1

## ----m2, eval=TRUE, echo=T, warning=FALSE, message=FALSE---------------------------------------------
res2<-blmpSDPD(formula = logc ~ logp+logy, data = data1, W = W,
               index = c("state","year"),
               model = list("ols","sar","sdm","sem","sdem","slx"), 
               effect = "time")

## ----m3, eval=TRUE, echo=T, warning=FALSE, message=FALSE---------------------------------------------
res3<-blmpSDPD(formula = logc ~ logp+logy, data = data1, W = W,
               index = c("state","year"),
               model = list("sar","sdm","sem","sdem"), 
               effect = "twoways",
               prior = "beta")

## ----m4, eval=TRUE, echo=T, warning=FALSE, message=FALSE---------------------------------------------
res4<-blmpSDPD(formula = logc ~ logp+logy, data = data1, W = W,
               index = c("state","year"),
               model = list("sar","sdm","sem","sdem","slx"), 
               effect = "twoways",
               ldet = "mc", ## log-determinant calculated with mcmc procedure
               dynamic = TRUE,
               prior = "uniform")

## ----m5, eval=TRUE, echo=T, warning=FALSE, message=FALSE---------------------------------------------
d2 <- plm::pdata.frame(data1, index=c('state', 'year'))
d2$llogc<-plm::lag(d2$logc) ## add lagged variable
data2<-d2[which(!is.na(d2$llogc)),]
rownames(data2)<-1:nrow(data2)
kk<-which(colnames(data2)=="llogc"); kk

res5<-blmpSDPD(formula = logc ~ logp+logy, data = data2, W = W,
               index = c("state","year"),
               model = list("sar","sdm","sem","sdem"), 
               effect = "individual",
               ldet = "full",
               dynamic = TRUE,
               tlaginfo = list(ind=kk),
               prior = "beta")

## ----mod1, eval=TRUE, echo=T, warning=FALSE, message=FALSE-------------------------------------------
mod1<-SDPDm(formula = logc ~ logp+logy, data = data1, W = W,
            index = c("state","year"),
            model = "sar", 
            effect = "individual")
summary(mod1)

## ---- eval=TRUE, echo=T, warning=FALSE, message=FALSE------------------------------------------------
mod2<-SDPDm(formula = logc ~ logp+logy, data = data1, W = W,
            index = c("state","year"),
            model = "sar", 
            effect = "individual",
            dynamic = T,
            tlaginfo = list(ind = NULL, tl = T, stl = T))
summary(mod2)

## ---- eval=TRUE, echo=T, warning=FALSE, message=FALSE------------------------------------------------
mod3<-SDPDm(formula = logc ~ logp+logy, data = data1, W = W,
            index = c("state","year"),
            model = "sar", 
            effect = "individual",
            LYtrans = T,
            dynamic = T,
            tlaginfo = list(ind = NULL, tl = T, stl = T))
summary(mod3)

## ---- eval=TRUE, echo=T, warning=FALSE, message=FALSE------------------------------------------------
mod4<-SDPDm(formula = logc ~ logp+logy, data = data2, W = W,
            index = c("state","year"),
            model = "sar", 
            effect = "individual",
            LYtrans = T,
            dynamic = T,
            tlaginfo = list(ind = kk, tl = T, stl = F))
summary(mod4)

## ---- eval=TRUE, echo=T, warning=FALSE, message=FALSE------------------------------------------------
mod5<-SDPDm(formula = logc ~ logp+logy, data = data1, W = W,
            index = c("state","year"),
            model = "sdm", 
            effect = "twoways",
            LYtrans = T,
            dynamic = T,
            tlaginfo = list(ind = NULL, tl = T, stl = T))
summary(mod5)

## ---- eval=TRUE, echo=T, warning=FALSE, message=FALSE------------------------------------------------
imp  <- impactsSDPDm(mod5)
imp

