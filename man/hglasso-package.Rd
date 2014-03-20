\name{hglasso-package}
\alias{hglasso-package}
\docType{package}
\title{
Fit the hub graphical lasso and hub covariance graph
}
\description{
This package is called hglasso, for "hub graphical lasso".  It implements two methods:hub graphical lasso and hub covariance graph.  All are described in the paper "Learning graphical models with hubs", by Tan et al. (2014).

The main functions are as follows:
(1) hglasso
(2) hcov

The first function, hglasso, performs hub graphical lasso.  The second function, hcov, performs hub covariance graph estimation.    

}
\details{
\tabular{ll}{
Package: \tab hglasso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-19\cr
License: \tab GPL (>=2.0) \cr
LazyLoad: \tab yes \cr
}

The package includes the following functinos:

\tabular{ll}{
\code{\link{hglasso}}: \tab Performs hub graphical lasso\cr
\code{\link{hcov}}:\tab Performs hub covariance graph estimation \cr
\code{\link{HubNetwork}}:\tab Generates inverse covariance matrix or covariance matrix with hubs \cr
\code{\link{image.hglasso}}: \tab Creates image plot of the matrix V and Z\cr
\code{\link{plot.hglasso}}: \tab Creates a graphical representation of the esitmated matrix Theta \cr
\code{\link{summary.hglasso}}: \tab Provides summary for the matrix Theta, Z, and V  \cr
}
}
\author{
Kean Ming Tan

Maintainer: Kean Ming Tan <keanming@uw.edu>
}
\references{
Tan et al. (2014). Learning graphical models with hubs.  arXiv.org/pdf/1402.7349.pdf.
}
\keyword{ package }

\seealso{
\code{\link{hglasso}}
\code{\link{hcov}}
}

\examples{
##############################################
# Example from Figure 1 in the manuscript
# A toy example to illustrate the results from 
# Hub Graphical Lasso
##############################################
#library(mvtnorm)
#set.seed(1)
#n=100
#p=100

# A network with 4 hubs
#network<-HubNetwork(p,0.99,4,0.1)
#Theta <- network$Theta
#truehub <- network$hubcol
# The four hub nodes have indices 14, 42, 45, 78
#print(truehub)

# Generate data matrix x
#x <- rmvnorm(n,rep(0,p),solve(Theta))
#x <- scale(x)

# Run Hub Graphical Lasso to estimate the inverse covariance matrix
# res1<-hglasso(cov(x),0.3,0.3,1.5)

# print out a summary of the object hglasso
#summary(res1)
# we see that the estimated hub nodes have indices 14, 42, 45, 78
# We successfully recover the 4 hub nodes

# Plot the matrices V and Z 
#image(res1)
#dev.off()
# Plot a graphical representation of the estimated inverse
# covariance matrix --- conditional independence graph
#plot(res1,main="Conditional Independence Graph")

}

