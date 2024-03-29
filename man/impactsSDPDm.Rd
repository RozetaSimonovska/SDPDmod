% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impactsSDPDm.R
\name{impactsSDPDm}
\alias{impactsSDPDm}
\title{Impacts for 'SDPDm' objects}
\usage{
impactsSDPDm(res, NSIM = 200, sd = 12345)
}
\arguments{
\item{res}{an object of class 'SDPDm'}

\item{NSIM}{number of simulations to be performed, default = 200}

\item{sd}{starting seed, default = 12345}
}
\value{
An object of class 'impactsSDPDm'
}
\description{
Direct and indirect effects estimates
}
\details{
For spatial dynamic panel data model:
\deqn{y_{t} = \tau y_{t-1} + \rho W y_{t} + \eta W y_{t-1} + X_{t} \beta + W X_{t} \theta + \alpha + \mu + u_{t}}
Short term effects for k\emph{th} explanatory variable:
\deqn{(I - \rho W)^{-1}(\beta_{k} I_{n} + \theta_{k} W)}
Long term effects for k\emph{th} explanatory variable:
\deqn{((1-\tau)I_{n} - (\rho+\eta)W)^{-1}(\beta_{k} I_{n} + \theta_{k} W)}
The direct effect is the average of the diagonal elements, and
the indirect effect is the average of the row sums of the non-diagonal elements
of the matrix.
}
\seealso{
\code{\link{SDPDm}}
}
\author{
Rozeta Simonovska
}
