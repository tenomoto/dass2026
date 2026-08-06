print_varnames <- function(x, dx = 0.5, ffact = 10, suffix="") {
  if (suffix == "sprd") {
    varnames <- c(paste("q0",suffix), paste("q2",suffix), paste("q4",suffix))
  } else {
    varnames <- c(paste("q0",suffix), paste("q2",suffix), paste("q4",suffix), paste("h",suffix), paste(ffact, "x f", suffix))
  }
  nc <- length(varnames)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), c(0, 1))
  for (v in varnames) {
    text(mean(x), 0.5, v, adj = 0.5)
    x <- x + 2 * pi + dx
  }
}

plot_waves_row <- function(x, q0, q2, q4, f, n, dx = 0.5, ffact = 10, add = FALSE, ...) {
  imax <- length(x)
  h <- calc_h(q0, q2, q4)
  ymat <- cbind(q0, q2, q4, h, ffact * f)
  ylim <- range(ymat)
  nc <- ncol(ymat)
  if (!add) {
    plot.new()
    plot.window(c(0, nc * (2 * pi + dx) - dx), ylim)
  }
  for (j in 1:nc) {
    if (n - 1 > 0 | j < nc) {
      lines(x, ymat[, j], ...)
      segments(x[1], 0, x[imax], 0)
    }
    x <- x + 2 * pi + dx
  }
  if (!add) {
    text(-0.3, 0, n - 1, xpd = TRUE, adj = 1)
  }
}

plot_waves <- function(x, hist, thist = NULL) {
  oldpar <- par(mfcol = c(tmax + 2, 1), xpd = TRUE, mar = c(0, 1, 0, 0))
  print_varnames(x)
  for (n in tmax:1) {
    plot_waves_row(x, hist$q0[, n], hist$q2[, n], hist$q4[, n], hist$f[, n], n)
    if (!is.null(thist))
    plot_waves_row(x, thist$q0[, n], thist$q2[, n], thist$q4[, n], thist$f[, n], n,
                   add = TRUE, lty = 2, col = "red")
  }
  par(oldpar)
}
