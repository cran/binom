\name{binom.probit}
\alias{binom.probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binomial confidence intervals using the probit parameterization}
\description{
  Uses the probit parameterization on the observed
  proportion to construct confidence intervals.
}
\usage{
binom.probit(x, n, conf.level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of number of successes in the binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence
    interval.}
  \item{\dots}{ignored}
}
\details{
  For derivations see \emph{doc/binom.pdf}.
}
\value{
  A \code{data.frame} containing the observed
  proportions and  the lower and upper bounds of the confidence
  interval.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.bayes}},
  \code{\link{binom.probit}}, \code{\link{binom.logit}}, \code{\link{binom.coverage}}}
\examples{
binom.probit(x = 0:10, n = 10)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
