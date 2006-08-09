\name{binom.profile}
\alias{binom.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binomial confidence intervals using the profile likelihood}
\description{
  Uses the profile likelihood on the observed
  proportion to construct confidence intervals.
}
\usage{
binom.profile(x, n, conf.level = 0.95, maxsteps = 50,
              del = zmax/5, bayes = TRUE, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of number of successes in the binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence
interval.}
\item{maxsteps}{The maximum number of steps to take in the profiles.}
\item{del}{The size of the step to take}
\item{bayes}{logical; if \code{TRUE} use a Bayesian correction at the
  edges.}
\item{plot}{logical; if \code{TRUE} plot the profile with a
  \code{spline} fit.}
\item{\dots}{ignored}
}
\details{
Confidence intervals are based on profiling the binomial deviance in the
neighbourhood of the MLE. If \code{x == 0} or \code{x == n} and
\code{bayes} is \code{TRUE}, then a Bayesian adjustment is made to move
the log-likelihood function away from \code{Inf}. Specifically, these
values are replaced by \code{(x + 0.5)/(n + 1)}, which is the posterier
mode of \code{f(p|x)} using Jeffrey's prior on \code{p}. Typically, the
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
}
\examples{
binom.profile(x = 0:10, n = 10)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
