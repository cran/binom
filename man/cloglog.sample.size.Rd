\name{cloglog.sample.size}
\alias{cloglog.sample.size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power and sample size}
\description{
Power and sample size for a binomial proportion using the cloglog parameterization. 
}
\usage{
cloglog.sample.size(p.alt, n = NULL, p = 0.5, power = 0.8, alpha = 0.05,
                    alternative = c("two.sided", "greater", "less"), exact.n = FALSE,
                    recompute.power = FALSE, phi = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.alt}{The alternative proportion in a one-sample test.}
  \item{n}{The sample size in a one-sample test. }
  \item{p}{The null proportion in a one-sample test. Default is 0.5. }
  \item{power}{The desired power level. Default is 0.80.}
  \item{alpha}{The desired alpha level - probability of a Type I
    error. Default is 0.05.
  }
  \item{alternative}{Nature of alternative hypothesis. One of
    "two.sided", "greater", "less".
  }
  \item{exact.n}{logical; If \code{TRUE}, the computed sample size will
    not be rounded up. Default is \code{FALSE}.
  }
  \item{recompute.power}{logical; If \code{TRUE}, after the sample size
    is computed, the power will be recomputed. This is only advantageous
    when the sample size is rounded up. Default is \code{FALSE}.
  }
  \item{phi}{Dispersion parameter by which to inflate (\code{phi > 1}) or
    deflate (\code{phi < 1}) variance. Default is 1.
  }
}
\details{
  This function can be used to calculate sample size, power or minimum
  detectable difference. It determines what to compute base on the
  arguments provided. If \code{p.alt} is given, but \code{n} is not, then
  sample size is computed. If \code{p.alt} is given along with \code{n},
  then the power is computed. If only \code{n} is provided, the minimum
  detectable difference is computed using the default power of 0.80.
}
\value{
  A \code{data.frame} containing the power, sample size and all of the
  input which was used to perform the computations.
}
\author{Sundar Dorai-Raj (sdorairaj@gmail.com)}
\seealso{\code{\link{binom.confint}}}
\examples{
cloglog.sample.size(p.alt = 0.8)
cloglog.sample.size(n = 20)
cloglog.sample.size(n = 20, power = 0.9)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{models}% __ONLY ONE__ keyword per line
