\name{binom.power}
\alias{binom.power}
\title{Power curves for binomial parameterizations}
\description{
  Uses Wald statistics to compute power curves for several parameterizations.
}
\usage{
binom.power(p.alt,
            n = 100,
            p = 0.5,
            alpha = 0.05,
            phi = 1,
            alternative = c("two.sided", "greater", "less"),
            method = c("cloglog", "logit", "probit", "asymp", "lrt", "exact"))
}
\arguments{
  \item{p.alt}{A vector of success probabilities under the alternative
    hypothesis.}
  \item{n}{A vector representing the number of independent trials in the
    binomial experiment.}
  \item{p}{A vector of success probabilities under the null hypothesis.}
  \item{alpha}{A vector of type-I error rates.}
  \item{phi}{A vector determining the overdispersion parameter for each
    binomial experiment.}
  \item{alternative}{Type of alternative hypothesis.}
  \item{method}{The method used to compute power.}
}
\details{
  For derivations see \emph{doc/binom.pdf}. \code{p.alt}, \code{n},
  \code{p}, \code{alpha}, and \code{phi} can all be vectors. The length
  of each argument will be expanded to the longest length. The function
  assumes the lengths are equal or can be wrapped for multiple values.
}
\value{
  The estimated probability of detecting the difference between
  \code{p.alt} and \code{p}.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.confint}}, \code{\link{binom.bayes}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}, \code{\link{binom.coverage}}}
\examples{
binom.power(0.95, alternative = "greater")
}
\keyword{univar}
\keyword{htest}
\keyword{models}
