binom.bayes <- function(x, n, conf.level = 0.95, type = c("highest", "central"),
                        prior.shape1 = 0.5, prior.shape2 = 0.5, # Jeffrey's prior
                        tol = .Machine$double.eps^.5, maxit = 1000, ...) {
  if(prior.shape1 <= 0 || prior.shape2 <= 0)
    stop("priors must be strictly positive.")
  if((length(x) != length(n)) && length(x) == 1) 
    x <- rep(x, length(n))
  if((length(x) != length(n)) && length(n) == 1) 
    n <- rep(n, length(x))
  ends <- x == 0 | x == n
  alpha <- rep(1 - conf.level, length(x))
  alpha[!ends] <- alpha[!ends] * 0.5
  a <- x + prior.shape1
  b <- n - x + prior.shape2
  p <- a/(a + b)
  lcl <- qbeta(alpha, a, b)
  ucl <- qbeta(1 - alpha, a, b)
  hpd <- any(pmatch(type, "highest", nomatch = FALSE))
  if(any(!ends) && hpd) {
    ci <- .C("binom_bayes",
             as.integer(x[!ends]),
             as.integer(n[!ends]),
             as.double(a[!ends]),
             as.double(b[!ends]),
             as.double(alpha[!ends]),
             lcl = as.double(lcl[!ends]),
             ucl = as.double(ucl[!ends]),
             as.integer(sum(!ends)),
             maxit = as.integer(maxit),
             tol = as.double(tol),
             PACKAGE = "binom")
    lcl[!ends] <- ci[["lcl"]]
    ucl[!ends] <- ci[["ucl"]]
  }
  lcl[x == 0] <- 0
  ucl[x == n] <- 1
  sig <- pbeta(lcl, a, b) + 1 - pbeta(ucl, a, b)
  res <- data.frame(method = "bayes", x = x, n = n, shape1 = a, shape2 = b,
                    mean = p, lower = lcl, upper = ucl, sig = sig)
  attr(res, "conf.level") <- conf.level
  res
}
