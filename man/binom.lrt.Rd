\name{binom.lrt}
\alias{binom.lrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binomial confidence intervals using the lrt likelihood}
\description{
  Uses the lrt likelihood on the observed
  proportion to construct confidence intervals.
}
\usage{
binom.lrt(x, n, conf.level = 0.95, bayes = FALSE, conf.adj = FALSE, plot
= FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of number of successes in the binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence
interval.}
\item{bayes}{logical; if \code{TRUE} use a Bayesian correction at the
  edges. Specfically, a beta prior with shape parameters 0.5 is used.
  If \code{bayes} is numeric, it is assumed to be the parameters
  to beta distribution.}
\item{conf.adj}{logical; if \code{TRUE} 0 or 100\% successes return a
  one-sided confidence interval}
\item{plot}{logical; if \code{TRUE} a plot showing the the square root
  of the binomial deviance with reference lines for mean, lower, and
  upper bounds. This argument can also be a \code{\link{list}} of plotting parameters
  to be passed to \code{\link[lattice]{xyplot}}.}
\item{\dots}{ignored}
}
\details{
Confidence intervals are based on profiling the binomial deviance in the
neighbourhood of the MLE. If \code{x == 0} or \code{x == n} and
\code{bayes} is \code{TRUE}, then a Bayesian adjustment is made to move
the log-likelihood function away from \code{Inf}. Specifically, these
values are replaced by \code{(x + 0.5)/(n + 1)}, which is the posterier
mode of \code{f(p|x)} using Jeffrey's prior on \code{p}. Furthermore, if
\code{conf.adj} is \code{TRUE}, then the upper (or lower) bound uses
a \code{1 - alpha} confidence level. Typically, the
observed mean will not be inside the estimated confidence interval. 
If \code{bayes} is \code{FALSE}, then the Clopper-Pearson exact method
is used on the endpoints. This tends to make confidence intervals at the
end too conservative, though the observed mean is guaranteed to be
within the estimated confidence limits.
}
\value{
  A \code{data.frame} containing the observed
  proportions and  the lower and upper bounds of the confidence
  interval.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.bayes}}, \code{\link{binom.cloglog}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}, \code{\link{binom.coverage}},
\code{\link[MASS]{confint}} in package MASS,
\code{\link[stats]{family}}, \code{\link[stats]{glm}}}
\examples{
binom.lrt(x = 0:10, n = 10)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
