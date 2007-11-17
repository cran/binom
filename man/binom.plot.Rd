\name{binom.plot}
\alias{binom.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage plots for binomial confidence intervals}
\description{
  Constructs coverage plots for binomial confidence intervals.
}
\usage{
binom.plot(n, method = binom.lrt, np = 500,
           conf.level = 0.95, actual = conf.level,
           type = c("xyplot", "levelplot"),
           tol = .Machine$double.eps^0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence
    interval.}
  \item{np}{Number of points to use in the plot.}
  \item{method}{The method used to estimate the confidence interval.}
  \item{actual}{The actual confidence interval used in the confidence
    interval. See Details.}
  \item{type}{See Details.}
  \item{tol}{The minimum probability of success to use in the plot.}
  \item{\dots}{Additional arguments to pass to
    \code{\link[lattice]{panel.xyplot}} or \code{\link[lattice]{panel.levelplot}}.}
}
\details{
  If \code{type} is \code{"xyplot"}, a line plot is created with
  coverage on the y-axis and binomial probability on the x-axis.
  A separate panel for every \code{n} is provided. If \code{actual} is
  provided then a horizontal reference line is added to the plot. This
  is only useful when \code{actual} is different from \code{conf.level},
  as is the case when calling \code{\link{binom.optim}}.

  If \code{type} is \code{"levelplot"}, a image plot is created with
  \code{x} = \code{0} to \code{n} on the vertical axis and binomial
  probability on the horizontal axis. Each row in the plot will be the
  confidence level for a given \code{x}. The color of the confidence
  interval is determined by the coverage probability. The argument
  \code{n} must only be of length one. If not, only the first \code{n}
  will be used and a warning is issued.

  In either plot type, the number of points at which the coverage
  probability is determined is specified by \code{np}. Increasing
  \code{np} gives a finer granularity but performance will suffer.
}
\value{
  An object of class \code{trellis}.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.optim}},
  \code{\link[lattice]{xyplot}}, \code{\link[lattice]{levelplot}}}
\examples{
binom.plot(5, type = "levelplot")
binom.plot(c(3, 5, 10, 25), type = "xyplot")
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{hplot}
