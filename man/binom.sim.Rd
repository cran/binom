\name{binom.sim}
\alias{binom.sim}
\title{Simulates confidence intervals for binomial data}
\description{
  Simulates binomial data for testing confidence interval coverage.
}
\usage{
binom.sim(M = 200, n = 100, p = 0.5, conf.level = 0.95, methods = "all",  ...)
}
\arguments{
  \item{M}{Number of simulations to create.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{p}{Probability of success under the null hypothesis.}
  \item{conf.level}{The level used in computing the confidence interval.}
  \item{methods}{The method used to compute power.}
  \item{\dots}{Additional arguments to pass to \code{\link{binom.confint}}}
}
\details{
  \code{M} binomial observations are created using \code{rbinom(M, n,
    p)}. The average number of times a confidence interval covers
  \code{p} is returned.
}
\value{
  The estimated coverage based on which method is requested.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.bayes}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}, \code{\link{binom.coverage}}}
\examples{
binom.sim(p = 0.95, method = "lrt")
}
\keyword{univar}
\keyword{htest}
\keyword{models}
