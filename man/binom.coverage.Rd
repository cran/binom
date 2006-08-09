\name{binom.coverage}
\alias{binom.coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probability coverage for binomial confidence intervals}
\description{
Determines the probability coverage for a binomial confidence interval. 
}
\usage{
binom.coverage(p, n, conf.level = 0.95, method = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{The (true) probability of success in a binomial experiment.}
  \item{n}{Vector of number of independent trials in the binomial experiment. }
  \item{conf.level}{The level of confidence to be used in the confidence interval. }
  \item{method}{Either a character string to be passed to
    \code{\link{binom.confint}} or a function that computes the upper
    and lower confidence bound for a binomial proportion. If a function
    is supplied, the first three arguments must be the same as
    \code{\link{binom.confint}} and the return value of the function
    must be a \code{data.frame} with column headers \code{"method"},
    \code{"lower"}, and \code{"upper"}. See \code{\link{binom.confint}}
    for available methods. Default is \code{"all"}.
  }
  \item{\dots}{Additional parameters to be passed to
    \code{\link{binom.confint}}. Only used when method is either
    \code{"bayes"} or \code{"profile"}
  }
}
\details{
  Derivations are based on the results given in the references. Methods
  whose coverage probabilities are consistently closer to 0.95 are more
  desireable. Thus, Wilson's, logit, and cloglog appear to be good for
  this sample size, while Jeffreys, asymptotic, and prop.test are
  poor. Jeffreys is a variation of Bayes using prior shape parameters of
  0.5 and having equal probabilities in the tail. The Jeffreys'
  equal-tailed interval was created using binom.bayes using (0.5,0.5) as
  the prior shape parameters and \code{type = "central"}.
}
\value{
  A \code{data.frame} containing the \code{"method"} used, \code{"n"}, \code{"p"},
  and the coverage probability, \code{C(p,n)}.
}
\references{
  L.D. Brown, T.T. Cai and A. DasGupta (2001), Interval
  estimation for a binomial proportion (with discussion), \emph{Statistical
  Science}, \bold{16}:101-133.
  
  L.D. Brown, T.T. Cai and A. DasGupta (2002), Confidence Intervals for
  a Binomial Proportion and Asymptotic Expansions, \emph{Annals of Statistics},
  \bold{30}:160-201.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com)}
\seealso{\code{\link{binom.confint}}, \code{\link{binom.length}}}
\examples{
binom.coverage(p = 0.5, n = 50)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
