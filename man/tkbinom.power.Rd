\name{tkbinom.power}
\alias{tkbinom.power}
\title{Power curves for binomial parameterizations}
\description{
  A Tcl/Tk graphics wrapper for \code{\link{binom.power}}.
}
\usage{
  tkbinom.power()
}
\details{
  A wrapper for \code{\link{binom.power}} that creates power curves
  based on user input.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }

\seealso{\code{\link{binom.power}}, \code{\link{binom.confint}}, \code{\link{binom.bayes}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}, \code{\link{binom.coverage}}}
\examples{
binom.power(0.95, alternative = "greater")
}
\keyword{univar}
\keyword{htest}
\keyword{models}
