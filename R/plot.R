plot_waves_overlay <- function(x, q0, q2, q4, f, i) {
  h <- calc_h(q0, q2, q4)
  qr <- range(c(q0, q2, q4, h, f))
  plot(x, q0, type = "l", ylim = c(qr[1], qr[2]), col = "blue",
       main = paste("step", i),
       xlab = "phase", ylab = "amplitude", xaxt = "n")
  axis(1, c(0, pi/2, pi, 3 * pi / 2, 2 * pi),
       c("0", "π/2", "π", "3π/2", "2π"))
  abline(h = 0, lty = 2)
  lines(x, q2, col = "red")
  lines(x, q4, col = "orange")
  lines(x, h, col = "black")
  lines(x, f, col = "purple")
  legend("topright", c("q0", "q2", "q4", "h", "f"), lty = 1,
         col = c("blue", "red", "orange", "black", "purple"))
}

print_varnames <- function(x, dx = 0.5, ffact = 10) {
  varnames <- c("q0", "q2", "q4", "h", paste(ffact, "x f"))
  nc <- length(varnames)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), c(0, 1))
  for (v in varnames) {
    text(mean(x), 0.5, v, adj = 0.5)
    x <- x + 2 * pi + dx
  }
}

plot_waves_row <- function(x, q0, q2, q4, f, n, dx = 0.5, ffact = 10) {
  h <- calc_h(q0, q2, q4)
  ymat <- cbind(q0, q2, q4, h, ffact * f)
  ylim <- range(ymat)
  nc <- ncol(ymat)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), ylim)
  for (j in 1:nc) {
    if (n - 1 > 0 | j < nc) {
      lines(x, ymat[, j])
      segments(x[1], 0, x[imax], 0)
    }
    x <- x + 2 * pi + dx
  }
  text(-0.3, 0, n - 1, xpd = TRUE, adj = 1)
}

plot_waves <- function(x, hist) {
  oldpar <- par(mfcol = c(tmax + 2, 1), xpd = TRUE, mar = c(0, 1, 0, 0))
  print_varnames(x)
  for (n in tmax:1) {
    plot_waves_row(x, hist$q0[, n], hist$q2[, n], hist$q4[, n], hist$f[, n], n)
  }
  par(oldpar)
}
