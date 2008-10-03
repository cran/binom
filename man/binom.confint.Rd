\name{binom.confint}
\alias{binom.agresti.coull}
\alias{binom.asymp}
\alias{binom.confint}
\alias{binom.exact}
\alias{binom.prop.test}
\alias{binom.wilson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binomial confidence intervals}
\description{
Uses eight different methods to obtain a confidence interval on the binomial probability. }
\usage{
binom.confint(x, n, conf.level = 0.95, methods = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of number of successes in the binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{conf.level}{The level of confidence to be used in the confidence interval.}
  \item{methods}{Which method to use to construct the interval. Any
    combination of \code{c("exact", "ac", "asymptotic", "wilson",
      "prop.test", "bayes", "logit", "cloglog", "probit")} is allowed. Default is
    \code{"all"}. }
  \item{\dots}{Additional arguments to be passed to \code{binom.bayes}.}
}
\details{
Nine methods are allowed for constructing the confidence interval(s):

\itemize{
\item{\code{exact}}{ - Pearson-Klopper method. See also
  \code{\link{binom.test}}.}
\item{\code{asymptotic}}{ - the text-book definition for confidence
  limits on a single proportion using the Central Limit Theorem.}
\item{\code{agresti-coull}}{ - Agresti-Coull method.}
\item{\code{wilson}}{ - Wilson method.}
\item{\code{prop.test}}{ - equivalent to \code{prop.test(x = x, n = n,
    conf.level = conf.level)$conf.int}.}
\item{\code{bayes}}{ - see \code{\link{binom.bayes}}.}
\item{\code{logit}}{ - see \code{\link{binom.logit}}.}
\item{\code{cloglog}}{ - see \code{\link{binom.cloglog}}.}
\item{\code{probit}}{ - see \code{\link{binom.probit}}.}
\item{\code{profile}}{ - see \code{\link{binom.profile}}.}
}

By default all eight are estimated for each value of \code{x} and/or
\code{n}. For the "logit", "cloglog", "probit", and "profile"
methods, the cases where \code{x == 0} or \code{x == n} are treated
separately. Specifically, the lower bound is replaced by
\code{(alpha/2)^n} and the upper bound is replaced by \code{(1-alpha/2)^n}.
}
\value{
  A \code{data.frame} containing the observed proportions and
  the lower and upper bounds of the confidence interval for all the
  methods in \code{"methods"}.
}
\references{
  A. Agresti and B.A. Coull (1998), Approximate is better than "exact"
  for interval estimation of binomial proportions, \emph{American
  Statistician}, \bold{52}:119-126.

  R.G. Newcombe, Logit confidence intervals and the inverse sinh
  transformation (2001), \emph{American Statistician}, \bold{55}:200-202.

  L.D. Brown, T.T. Cai and A. DasGupta (2001), Interval estimation for a
  binomial proportion (with discussion), \emph{Statistical Science},
  \bold{16}:101-133.

  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (1997)
  \emph{Bayesian Data Analysis}, London, U.K.: Chapman and Hall.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.bayes}}, \code{\link{binom.logit}},
  \code{\link{binom.probit}}, \code{\link{binom.cloglog}},
  \code{\link{binom.coverage}}, \code{\link[stats]{prop.test}},
  \code{\link[stats]{binom.test}} for comparison to method
  \code{"exact"}}
\examples{
binom.confint(x = c(2, 4), n = 100, tol = 1e-8)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
