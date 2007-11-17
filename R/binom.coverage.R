binom.coverage <- function(p, n, conf.level = 0.95, method = "all", ...) {
  if(missing(p)) p <- seq(0, 1, length = 200)
  x <- unlist(lapply(lapply(n, ":", 0), rev))
  n <- rep(n, n + 1)
  ci <- if(is.function(method)) {
    method(x, n, conf.level, ...)
  } else if(is.character(method) && exists(method) && method != "all") {
    get(method)(x, n, conf.level, ...)
  } else {
    binom.confint(x, n, conf.level, method, ...)
  }
  ci <- ci[c("method", "x", "n", "lower", "upper")]
  z <- merge(ci, data.frame(p = p))
  z$coverage <- with(z, (p >= lower & p <= upper) * dbinom(x, n, p))
  z <- aggregate(z["coverage"], z[c("method", "p", "n")], sum)
  z <- z[order(z$method, z$p), ]
  row.names(z) <- seq(NROW(z))
  z
}
