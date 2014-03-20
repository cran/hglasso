\name{summary.hglasso}
\alias{summary.hglasso}
\alias{summary.hcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot an object of class \code{hglasso} or \code{hcov}
}
\description{
This function provides some information for an object \code{\link{hglasso}} or \code{\link{hcov}}
}
\usage{
\method{summary}{hglasso}(object, \dots) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{\link{hglasso}} or \code{\link{hcov}}.
}

  \item{\dots}{
any other arguments passed to \code{print}.
}

}
\details{
Some information for an object hglasso or hcov: (1) The number of observations n and the number of features p.  (2) The number of edges in Theta, V, and Z. (3) The indices for hub nodes, and also the number of edges within each hub node.}

\references{
Tan et al. (2014). Learning graphical models with hubs.  arXiv.org/pdf/1402.7349.pdf.
}
\author{
Kean Ming Tan 
}


\seealso{
\code{\link{image.hglasso}}
\code{\link{plot.hglasso}}
\code{\link{hglasso}}
\code{\link{hcov}}
}
\examples{
# See example in hglasso or hcov.
}

