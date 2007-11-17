\name{binom.power}
\alias{binom.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power curves for binomial parameterizations}
\description{
  Uses Wald statistics to compute power curves for several parameterizations.
}
\usage{
binom.power(p.alt, n = 100, p = 0.5, alpha = 0.05, phi = 1,
            alternative = c("two.sided", "greater", "less"),
            method = c("cloglog", "logit", "probit", "asymp", "lrt", "exact")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.alt}{Probability of success under the alternative hypothesis.}
  \item{n}{Vector of number of independent trials in the binomial experiment.}
  \item{p}{Probability of success under the null hypothesis.}
  \item{alpha}{Type-I error rate.}
  \item{phi}{Overdispersion parameter.}
  \item{alternative}{Type of alternative hypothesis.}
  \item{method}{The method used to compute power.}
}
\details{
  For derivations see \emph{doc/binom.pdf}.
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
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
