\name{binom.bayes}
\alias{binom.bayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binomial confidence intervals using Bayesian inference}
\description{
  Uses a beta prior on the probability of success for a binomial
  distribution, determines a two-sided confidence interval from a beta
  posterior.
}
\usage{
binom.bayes(x, n, conf.level = 0.95, type = c("highest", "central"),
            prior.shape1 = 0.5, prior.shape2 = 0.5,
            tol = .Machine$double.eps^0.5, maxit = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of number of successes in the binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence interval.}
  \item{type}{The type of confidence interval (see Details).}
  \item{prior.shape1}{The value of the first shape parameter to be used in the prior beta.}
  \item{prior.shape2}{The value of the second shape parameter to be used in the prior beta.}
  \item{tol}{A tolerance to be used in determining the highest probability density interval.}
  \item{maxit}{Maximum number of iterations to be used in determining
    the highest probability interval.}
  \item{\dots}{Ignored.}
}
\details{
  Using the conjugate beta prior on the distribution of p (the
  probability of success) in a binomial experiment, constructs a
  confidence interval from the beta posterior. From Bayes theorem the
  posterior distribution of p given the data x is:
  
  \code{p|x ~ Beta(x + prior.shape1, n - x + prior.shape2)}

  The default prior is Jeffrey's prior which is a Beta(0.5, 0.5)
  distribution. Thus the posterior mean is \code{(x + 0.5)/(n + 1)}.
  
  The default type of interval constructed is "highest" which computes
  the highest probability density (hpd) interval which assures the
  shortest interval possible. The hpd intervals will achieve a
  probability that is within tol of the specified conf.level. Setting
  type to "central" constructs intervals that have equal tail
  probabilities.
  
  If 0 or n successes are observed, a one-sided confidence interval is
  returned.
}
\value{
  A data.frame containing the observed proportions and
  the lower and upper bounds of the confidence interval.
}
\references{
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (1997)
  \emph{Bayesian Data Analysis}, London, U.K.: Chapman and Hall.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }
\seealso{\code{\link{binom.confint}}, \code{\link{binom.cloglog}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}}
\examples{
binom.bayes(x = 0:10, n = 10, tol = 1e-9)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
