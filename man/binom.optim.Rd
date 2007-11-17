\name{binom.optim}
\alias{binom.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal binomial confidence intervals}
\description{
  Uses optimization to minimize the integrated mean squared error
  between the calculated coverage and the desired confidence level
  for a given binomial confidence interval.
}
\usage{
binom.optim(n, conf.level = 0.95, method = binom.lrt,
            k = n\%/\%2 + 1, p0 = 0, transform = TRUE,
            plot = FALSE, tol = .Machine$double.eps^0.5,
            start = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence
    interval.}
  \item{method}{The method used to estimate the confidence interval.}
  \item{k}{See Details.}
  \item{p0}{The minimum probability of success to allow in the
    optimization. See Details.}
  \item{transform}{logical; If \code{TRUE} the optimizer will do an
    unconstrained optimization on the signficance probability in the
    logit space.}
  \item{plot}{logical; If \code{TRUE} the results are sent to \code{\link{binom.plot}}.}
  \item{tol}{The minimum significance level to allow in the
    optimization. See Details.}
  \item{start}{A starting value on the optimal confidence level.}
  \item{\dots}{Additional arguments to pass to \code{\link[stats]{optim}}.}
}
\details{
This function minimizes the squared error between the expected coverage
probability and the desired confidence level.

\deqn{\alpha_{opt}=\arg\min_{\alpha}\int_{0}^{1}[C(p,n)-(1-\alpha)^2dp}{%
     alpha[opt]=argmin[alpha] integral[C(p,n)-(1-alpha)]^2dp}

The optimizer will adjust confidence intervals for all \code{x} =
\code{0} to \code{n} depending on the value of \code{k} provided. If
\code{k} is one, only the confidence levels for \code{x} = \code{0} and
\code{n} are adjusted. If \code{k} = \code{[n/2]} then all confidence
intervals are adjusted. This assumes the confidence intervals are the
same length for \code{x} = \code{x[k]} and \code{x[n - k + 1]}, which is
the case for all methods provided in this package except
\code{\link{binom.cloglog}}.}
\value{
  A \code{list} with the following elements:

  \itemize{
    \item{par}{Final confidence levels. The length of this vector is
      \code{k}.}
    \item{value}{The final minimized value from \code{\link{optim}}.}
    \item{counts}{The number of function and gradient calls from
      \code{\link{optim}}.}
    \item{convergence}{Convergence code from \code{\link{optim}}.}
    \item{message}{Any message returned by the L-BFGS-B or BFGS
      optimizer.}
    \item{confint}{A \code{data.frame} returned from a call to
      \code{method} using the optimized confidence levels.}
  }
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.plot}},
  \code{\link{binom.coverage}}, \code{\link[stats]{optim}}}
\examples{
binom.optim(10, k = 1) ## determine optimal significance for x = 0, 10 only
binom.optim(3, method = binom.wilson) ## determine optimal significance for all x
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{optimize}
