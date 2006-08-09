\name{tkbinom.power}
\alias{tkbinom.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power curves for binomial parameterizations}
\description{
  A Tcl/Tk graphics wrapper for \code{\link{binom.power}}.
}
\usage{
  tkbinom.power()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\details{
  A wrapper for \code{\link{binom.power}} that creates power curves
  based on user input.
}
\value{
  None.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.power}}, \code{\link{binom.confint}}, \code{\link{binom.bayes}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}, \code{\link{binom.coverage}}}
}
\examples{
binom.power(0.95, alternative = "greater")
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
