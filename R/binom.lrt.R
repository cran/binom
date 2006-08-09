binom.lrt <- function(x, n, conf.level = 0.95, bayes = FALSE, conf.adj = FALSE, plot = FALSE, ...) {
  do.plot <- (is.logical(plot) && plot) || is.list(plot)
  if(do.plot) require(lattice, TRUE)
  xn <- cbind(x = x, n = n)
  ok <- !is.na(xn[, 1]) & !is.na(xn[, 2])
  x <- xn[ok, "x"]
  n <- xn[ok, "n"]
  p <- ifelse(ok, x/n, NA)
  res <- data.frame(xn, mean = p)
  alpha <- 1 - conf.level
  alpha <- rep(alpha, length = length(p))
  res$lower <- rep(0, NROW(res))
  res$upper <- rep(1, NROW(res))
  args <- list(...)
  tol <- if(is.null(args$tol)) .Machine$double.eps^0.5 else args$tol
  bindev <- function(y, x, mu, wt, bound = 0, tol = .Machine$double.eps^0.5, ...) {
    ll.y <- ifelse(y %in% c(0, 1), 0, binom:::ldbinom(x, wt, y))
    ll.mu <- ifelse(mu %in% c(0, 1), 0, binom:::ldbinom(x, wt, mu))
    f <- ifelse(abs(y - mu) < tol, 0, sign(y - mu) * sqrt(-2 * (ll.y - ll.mu)))
    f - bound
  }
  args$f <- bindev
  x0 <- x == 0
  xn <- x == n
  z <- qnorm(1 - alpha * ifelse((x0 | xn) & conf.adj, 0.25, 0.5))
  if((is.logical(bayes) && bayes) || is.numeric(bayes)) {
    ## use bayesian adjustment for edge correction
    if(is.logical(bayes)) {
      bayes <- c(0.5, 0.5)
    } else if(length(bayes) == 1) {
      bayes <- c(bayes, bayes)
    }
    if(any(edge <- x0 | xn)) {
      n[edge] <- n[edge] + bayes[1] + bayes[2]
      x[edge] <- x[edge] + bayes[1]
    }
    p <- x/n
    bayes <- TRUE
  }
  plot.df <- list()
  for(i in seq(NROW(res))) {
    if(!ok[i]) {
      res$lower[i] <- res$upper[i] <- NA
      next
    }
    args[c("x", "mu", "wt")] <- list(x = x[i], mu = p[i], wt = n[i])
    if(!x0[i] && tol < p[i]) {
      args$interval <- c(tol, if(p[i] < tol || p[i] == 1) 1 - tol else p[i])
      args$bound <- -z[i]
      res$lower[i] <- if(bindev(tol, x[i], p[i], n[i], -z[i], tol) > 0) {
        if(conf.adj) z[i] <- qnorm(1 - alpha)
        0
      } else {
        do.call("uniroot", args)$root
      }
    }
    if(!xn[i] && p[i] < 1 - tol) {
      args$interval <- c(if(p[i] > 1 - tol) tol else p[i], 1 - tol)
      args$bound <- z[i]
      res$upper[i] <- if(bindev(1 - tol, x[i], args$interval[1], n[i], z[i], tol) < 0) {
        args$interval <- c(tol, p[i])
        if(conf.adj) z[i] <- qnorm(1 - alpha)
        args$bound <- -z[i]
        res$lower[i] <- if(bindev(tol, x[i], if(p[i] < tol || p[i] == 1) 1 - tol else p[i], n[i], -z[i], tol) > 0) {
          0
        } else {
          do.call("uniroot", args)$root
        }
        1
      } else {
        do.call("uniroot", args)$root
      }
    }
    if(do.plot) {
      min.p <- res$lower[i]
      max.p <- res$upper[i]
      dx.p <- 0.1 * (max.p - min.p)
      expand.p <- seq(max(min.p - dx.p, tol), min(max.p + dx.p, 1 - tol), len = 100)
      plot.df[[i]] <- data.frame(p = expand.p)
      plot.df[[i]]$x <- res$x[i]
      plot.df[[i]]$n <- res$n[i]
      plot.df[[i]]$lower <- res$lower[i]
      plot.df[[i]]$upper <- res$upper[i]
      plot.df[[i]]$mean <- res$mean[i]
      plot.df[[i]]$deviance <- bindev(expand.p, x[i], p[i], n[i], 0, tol)
    }
  }
  if(do.plot) {
    plot.df <- do.call("rbind", plot.df)
    xn <- unique(data.frame(x = plot.df$x, n = plot.df$n))
    xn <- xn[order(xn$x, xn$n), ]
    lev.xn <- with(xn, sprintf("x = %s; n = %s", x, n))
    plot.df$xn <- with(plot.df, factor(sprintf("x = %s; n = %s", x, n), levels = lev.xn))
    plot.args <- if(is.list(plot)) plot else list()
    plot.args$x <- deviance ~ p | xn
    plot.args$data <- plot.df
    plot.args$lower <- plot.df$lower
    plot.args$upper <- plot.df$upper
    plot.args$mean <- plot.df$mean
    plot.args$alpha <- unique(alpha)
    if(is.null(plot.args$panel))
      plot.args$panel <- panel.binom.lrt
#    if(is.null(plot.args$scales)) {
#      plot.args$scales <- list(x = list(relation = "free"))
#    } else if(is.null(args$plot.args$scales$x)) {
#      plot.args$scales$x <- list(relation = "free")
#    }
    if(is.null(plot.args$xlab))
      plot.args$xlab <- "Probability of Success (p)"
    if(is.null(plot.args$ylab))
      plot.args$ylab <- "Standard Normal Quantiles"
    if(is.null(plot.args$par.strip.text)) {
      plot.args$par.strip.text <- list(cex = 1.2)
    } else if(is.null(plot.args$par.strip.text)) {
      plot.args$par.strip.text$cex <- 1.2
    }
    if(is.null(plot.args$as.table)) plot.args$as.table <- TRUE
    print(do.call("xyplot", plot.args))
  }
  attr(res, "conf.level") <- conf.level
  cbind(method = "lrt", res)
}

panel.binom.lrt <- function(x, y, subscripts, lower, upper, mean, alpha, ...) {
  s <- spline(x, y)
  panel.xyplot(s$x, s$y, type = "l", lwd = 2, col = "darkblue")
  lower <- unique(lower[subscripts])
  upper <- unique(upper[subscripts])
  mean <- unique(mean[subscripts])
  limits <- current.panel.limits()
  xlim <- limits$xlim
  ylim <- limits$ylim
  adj.x <- diff(xlim) * 0.04
  adj.y <- diff(ylim) * 0.04
  llines(c(xlim[1], mean), c(0, 0), lty = 2, col = "darkgray", lwd = 2)
  llines(c(mean, mean), c(ylim[1], 0), lty = 2, col = "darkgray", lwd = 2)
  lab <- sprintf("p = %0.4f", mean)
  ltext(xlim[1] + adj.x, adj.y, lab, adj = 0, col = "black")
  lpoints(mean, 0, col = "darkred", pch = 16, cex = 1.2)
  if(lower > xlim[1] && lower < mean) {
    q <- qnorm(alpha/2)
    llines(c(xlim[1], lower), c(q, q), lty = 2, col = "darkgray", lwd = 2)
    llines(c(lower, lower), c(ylim[1], q), lty = 2, col = "darkgray", lwd = 2)
    ##lab <- substitute(hat(p)[L] == phat, list(phat = sprintf("%0.4f", lower)))
    lab <- sprintf("LCL = %0.4f", lower)
    ltext(xlim[1] + adj.x, q + adj.y, lab, adj = 0, col = "black")
    lpoints(lower, q, col = "darkred", pch = 16, cex = 1.2)
  }
  if(upper < xlim[2] && mean < upper) {
    q <- qnorm(1 - alpha/2)
    llines(c(xlim[1], upper), c(q, q), lty = 2, col = "darkgray", lwd = 2)
    llines(c(upper, upper), c(ylim[1], q), lty = 2, col = "darkgray", lwd = 2)
    ##lab <- substitute(hat(p)[U] == phat, list(phat = sprintf("%0.4f", upper)))
    lab <- sprintf("UCL = %0.4f", upper)
    ltext(xlim[1] + adj.x, q + adj.y, lab, adj = 0, col = "black")
    lpoints(upper, q, col = "darkred", pch = 16, cex = 1.2)
  }
}


#binom.loglik <- function(x, n, conf.level = 0.95, bayes = TRUE, conf.adj = FALSE, edge.adj = TRUE, ...) {
#  if(length(x) > 1 || length(n) > 1)
#    stop("binom.loglik is not vectorised")
#  ok <- !is.na(x) & !is.na(n) & x >=0 & x <= n & n > 0
#  if(!ok) return(c(lower = NA, upper = NA))
#  p <- x/n
#  alpha <- 1 - conf.level
#  lower <- 0
#  upper <- 1
#  args <- list(...)
#  tol <- if(is.null(args$tol)) .Machine$double.eps else args$tol
#  args$f <- function(y, mu, w, bound) {
##    x <- p * n
##    lr0 <- ifelse(x == 0, 0, log(p0/p))
##    lr1 <- ifelse(x == n, 0, log((1 - p0)/(1 - p)))
##    lr <- -2 * (x * lr0 + (n - x) * lr1)
##    sign(p0 - p) * sqrt(lr) - bound
#    devy <- y
#    nz <- y != 0
#    devy[nz] <- y[nz] * log(y[nz])
#    nz <- (1 - y) != 0
#    devy[nz] <- devy[nz] + (1 - y[nz]) * log(1 - y[nz])
#    devmu <- y * log(mu) + (1 - y) * log(1 - mu)
#    if(any(small <- mu * (1 - mu) < .Machine$double.eps)) {
#      smu <- mu[small]
#      sy <- y[small]
#      smu <- ifelse(smu < .Machine$double.eps, .Machine$double.eps, smu)
#      onemsmu <- ifelse((1 - smu) < .Machine$double.eps, .Machine$double.eps, 1 - smu)
#      devmu[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
#    }
#    devi <- 2 * (devy - devmu)
#    sign(y - mu) * sqrt(sum(w * devi)) - bound
#  }
#  x0 <- x < 1
#  xn <- x > n - 1
#  z <- qnorm(1 - alpha * if((x0 || xn) && conf.adj) 1 else 0.5)
#  if((is.logical(bayes) && bayes) || is.numeric(bayes)) {
#    if(is.logical(bayes)) {
#      bayes <- c(0.5, 0.5)
#    } else if(length(bayes) == 1) {
#      bayes <- c(bayes, bayes)
#    }
#    if(x0 || xn) {
#      n <- n + bayes[1] + bayes[2]
#      x <- x + bayes[1]
#    }
#    p <- x/n
#    bayes <- TRUE
#  }
#  if((x0 || xn) && !bayes) {
#    if(xn) {
#      lower <- (alpha/2)^(1/n)
#    } else {
#      upper <- 1 - (alpha/2)^(1/n)
#    }
#  } else {
#    args[c("mu", "w")] <- list(mu = p, w = n)
#    if(tol > p || p > 1 - tol)
#      warning("estimated p is less than tolerance for uniroot")
#    if((!x0 || edge.adj) && tol < p) {
#      args$interval <- c(tol, p)
#      args$bound <- -z
#      lower <- if(args$f(tol, p, n, -z) > 0) 0 else do.call("uniroot", args)$root
#    }
#    if((!xn || edge.adj) && p < 1 - tol) {
#      args$interval <- c(p, 1 - tol)
#      args$bound <- z
#      upper <- if(args$f(1 - tol, p, n, z) < 0) 1 else do.call("uniroot", args)$root
#    }
#  }
#  c(lower = lower, upper = upper)
#}
