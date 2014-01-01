\name{binom.bayes}
\alias{binom.bayes}
\alias{binom.bayes.densityplot}
\title{Binomial confidence intervals using Bayesian inference}
\description{
  Uses a beta prior on the probability of success for a binomial
  distribution, determines a two-sided confidence interval from a beta
  posterior. A plotting function is also provided to show the
  probability regions defined by each confidence interval.
}
\usage{
binom.bayes(x, n,
            conf.level = 0.95,
            type = c("highest", "central"),
            prior.shape1 = 0.5,
            prior.shape2 = 0.5,
            tol = .Machine$double.eps^0.5,
            maxit = 1000, ...)

binom.bayes.densityplot(bayes,
                        npoints = 500,
                        fill.central = "lightgray",
                        fill.lower = "steelblue",
                        fill.upper = fill.lower,
                        alpha = 0.8, ...)
}
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
  \item{bayes}{The output \code{data.frame} from \code{binom.bayes}.}
  \item{npoints}{The number of points to use to draw the density
    curves. Higher numbers give smoother densities.}
  \item{fill.central}{The color for the central density.}
  \item{fill.lower,fill.upper}{The color(s) for the upper and lower
    density.}
  \item{alpha}{The alpha value for controlling transparency.}
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
  For \code{binom.bayes}, a \code{data.frame} containing the observed
  proportions and the lower and upper bounds of the confidence interval.

  For \code{binom.bayes.densityplot}, a \code{ggplot} object that can
  printed to a graphics device, or have additional layers added.
}
\references{
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (1997)
  \emph{Bayesian Data Analysis}, London, U.K.: Chapman and Hall.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com) }
\seealso{\code{\link{binom.confint}}, \code{\link{binom.cloglog}},
  \code{\link{binom.logit}}, \code{\link{binom.probit}}}
\examples{
# Example using highest probability density.
hpd <- binom.bayes(
  x = 0:10, n = 10, type = "highest", conf.level = 0.8, tol = 1e-9)
print(hpd)
binom.bayes.densityplot(hpd)
# Remove the extremes from the plot since they make things hard
# to see.
binom.bayes.densityplot(hpd[hpd$x != 0 & hpd$x != 10, ])

# Example using central probability.
central <- binom.bayes(
  x = 0:10, n = 10, type = "central", conf.level = 0.8, tol = 1e-9)
print(central)
binom.bayes.densityplot(central)
# Remove the extremes from the plot since they make things hard
# to see.
binom.bayes.densityplot(central[central$x != 0 & central$x != 10, ])

}
\keyword{univar}
\keyword{htest}
\keyword{models}
