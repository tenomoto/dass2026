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

plot_waves_with_truth_row <- function(x, q0, q2, q4, f, q0_t, q2_t, q4_t, f_t, n, dx = 0.5, ffact = 10) {
  h <- calc_h(q0, q2, q4)
  ymat <- cbind(q0, q2, q4, h, ffact * f)
  h_t <- calc_h(q0_t, q2_t, q4_t)
  ymat_t <- cbind(q0_t, q2_t, q4_t, h_t, ffact * f_t)
  ylim <- range(ymat)
  nc <- ncol(ymat)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), ylim)
  for (j in 1:nc) {
    if (n - 1 > 0 | j < nc) {
      lines(x, ymat[, j])
      lines(x, ymat_t[, j], col='red', lty=3)
      segments(x[1], 0, x[imax], 0)
    }
    x <- x + 2 * pi + dx
  }
  text(-0.3, 0, n - 1, xpd = TRUE, adj = 1)
}

plot_waves_with_truth <- function(x, hist, truth, suffix="") {
  oldpar <- par(mfcol = c(tmax + 2, 1), xpd = TRUE, mar = c(0, 1, 0, 0))
  print_varnames(x,suffix=suffix)
  for (n in tmax:1) {
    plot_waves_with_truth_row(x, 
      hist$q0[, n], hist$q2[, n], hist$q4[, n], hist$f[, n], 
      truth$q0[, n], truth$q2[, n], truth$q4[, n], truth$f[, n], 
      n)
  }
  par(oldpar)
}

plot_diff_row <- function(x, q0, q2, q4, f, q0_t, q2_t, q4_t, f_t, n, dx = 0.5, ffact = 10, ylim = NA) {
  h <- calc_h(q0, q2, q4)
  h_t <- calc_h(q0_t, q2_t, q4_t)
  ymat <- cbind(q0-q0_t, q2-q2_t, q4-q4_t, h-h_t, ffact * (f-f_t))
  if (is.na(ylim[1])|is.na(ylim[2])) {
    ylim <- range(ymat)
  }
  nc <- ncol(ymat)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), ylim)
  for (j in 1:nc) {
    if (n - 1 > 0 | j < nc) {
      lines(x, ymat[, j], col='red')
      segments(x[1], 0, x[imax], 0)
    }
    x <- x + 2 * pi + dx
  }
  text(-0.3, 0, n - 1, xpd = TRUE, adj = 1)
}

plot_diff <- function(x, hist, truth, ffact=10.0) {
  oldpar <- par(mfcol = c(tmax + 2, 1), xpd = TRUE, mar = c(0, 1, 0, 0))
  print_varnames(x,suffix="error",ffact=ffact)
  q0 <- hist$q0[, 1]; q2 <- hist$q2[, 1]; q4 <- hist$q4[, 1]
  q0_t <- truth$q0[, 1]; q2_t <- truth$q2[, 1]; q4_t <- truth$q4[, 1]
  h <- calc_h(q0, q2, q4)
  h_t <- calc_h(q0_t, q2_t, q4_t)
  f <- hist$f[, 2]
  f_t <- truth$f[, 2]
  ymat <- cbind(q0-q0_t, q2-q2_t, q4-q4_t, h-h_t, ffact*(f-f_t))
  ylim <- range(ymat)
  for (n in tmax:1) {
    plot_diff_row(x, 
                  hist$q0[, n], hist$q2[, n], hist$q4[, n], hist$f[, n], 
                  truth$q0[, n], truth$q2[, n], truth$q4[, n], truth$f[, n], 
                  n,ylim=ylim)
  }
  par(oldpar)
}

plot_sprd_row <- function(x, q0s, q2s, q4s, n, dx = 0.5, ffact = 10, ylim = c(NA,NA)) {
  ymat <- cbind(q0s, q2s, q4s)
  if (is.na(ylim[1])|is.na(ylim[2])) {
    ylim <- range(ymat)
    ylim[1] <- 0.0
  }
  nc <- ncol(ymat)
  plot.new()
  plot.window(c(0, nc * (2 * pi + dx) - dx), ylim)
  for (j in 1:nc) {
    lines(x, ymat[, j], col='red')
    segments(x[1], 0, x[imax], 0)
    x <- x + 2 * pi + dx
  }
  text(-0.3, 0, n - 1, xpd = TRUE, adj = c(1,0))
}

plot_sprd <- function(x, hist, dx = 0.5) {
  oldpar <- par(mfcol = c(tmax + 2, 1), xpd = TRUE, mar = c(0, 1, 0, 0))
  print_varnames(x, suffix="sprd")
  ymat <- cbind(hist$q0s[, 1], hist$q2s[, 1], hist$q4s[, 1])
  ylim <- range(ymat)
  ylim[1] <- 0.0
  for (n in tmax:1) {
    plot_sprd_row(x, hist$q0s[, n], hist$q2s[, n], hist$q4s[, n], n, ylim=ylim)
  }
  par(oldpar)
}

calc_rmse <- function(hist, truth, guess) {
  h <- calc_h(hist$q0, hist$q2, hist$q4)
  h_t <- calc_h(truth$q0, truth$q2, truth$q4)
  h_g <- calc_h(guess$q0, guess$q2, guess$q4)
  rmse_q0 <- sqrt(colMeans((hist$q0 - truth$q0)^2))
  rmse_q2 <- sqrt(colMeans((hist$q2 - truth$q2)^2))
  rmse_q4 <- sqrt(colMeans((hist$q4 - truth$q4)^2))
  rmse_h <- sqrt(colMeans((h - h_t)^2))
  grmse_q0 <- sqrt(colMeans((guess$q0 - truth$q0)^2))
  grmse_q2 <- sqrt(colMeans((guess$q2 - truth$q2)^2))
  grmse_q4 <- sqrt(colMeans((guess$q4 - truth$q4)^2))
  grmse_h <- sqrt(colMeans((h_g - h_t)^2))
  rmse <- data.frame(q0 = rmse_q0, q2 = rmse_q2, q4 = rmse_q4, h = rmse_h,
               gq0 = grmse_q0, gq2 = grmse_q2, gq4 = grmse_q4, gh = grmse_h)
  if ("q0s" %in% names(hist)) {
    sprd_q0 <- colMeans(hist$q0s)
    sprd_q2 <- colMeans(hist$q2s)
    sprd_q4 <- colMeans(hist$q4s)
    rmse$q0s <- sprd_q0
    rmse$q2s <- sprd_q2
    rmse$q4s <- sprd_q4
  }
  rmse
}

plot_rmse <- function(rmse) {
  rmse_q0 <- rmse$q0
  rmse_q2 <- rmse$q2
  rmse_q4 <- rmse$q4
  rmse_h <- rmse$h
  grmse_q0 <- rmse$gq0
  grmse_q2 <- rmse$gq2
  grmse_q4 <- rmse$gq4
  grmse_h <- rmse$gh
  ylim <- range(c(rmse_q0,rmse_q2,rmse_q4,rmse_h,grmse_q0,grmse_q2,grmse_q4,grmse_h))
  if ("q0s" %in% names(rmse)) {
    sprd_q0 <- rmse$q0s
    sprd_q2 <- rmse$q2s
    sprd_q4 <- rmse$q4s
    ylim2 <- range(c(sprd_q0,sprd_q2,sprd_q4))
    ylim[1] <- min(ylim[1],ylim2[1])
    ylim[2] <- max(ylim[2],ylim2[2])
  }
  oldpar <- par(mfcol=c(1,1), mar=c(3,2,1,1))
  plot(1:length(rmse_q0), rmse_q0, type="l", col="black", lwd=2, ylim=ylim,
       ylab='RMSE')
  lines(1:length(rmse_q2), rmse_q2, col="red", lwd=2, ylim=ylim)
  lines(1:length(rmse_q4), rmse_q4, col="blue", lwd=2, ylim=ylim)
  lines(1:length(rmse_h), rmse_h, col="green", lwd=2, ylim=ylim)
  lines(1:length(grmse_q0), grmse_q0, col="black", lty=2, ylim=ylim)
  lines(1:length(grmse_q2), grmse_q2, col="red", lty=2, ylim=ylim)
  lines(1:length(grmse_q4), grmse_q4, col="blue", lty=2, ylim=ylim)
  lines(1:length(grmse_h), grmse_h, col="green", lty=2, ylim=ylim)
  if ("q0s" %in% names(rmse)) {
    lines(1:length(sprd_q0), sprd_q0, col="black", lty=3, ylim=ylim)
    lines(1:length(sprd_q2), sprd_q2, col="red", lty=3, ylim=ylim)
    lines(1:length(sprd_q4), sprd_q4, col="blue", lty=3, ylim=ylim)
    legend("topright",
         c("q0 DA","q2 DA","q4 DA","h DA",
           "q0 guess","q2 guess","q4 guess","h guess",
           "q0 sprd","q2 sprd","q4 sprd"), 
         col = rep(c("black","red","blue","green"),3), 
         lty = c(1,1,1,1,2,2,2,2,3,3,3), 
         lwd = c(2,2,2,2,1,1,1,1,1,1,1), ncol=3)
  } else {
    legend("topright",
         c("q0 DA","q2 DA","q4 DA","h DA",
           "q0 guess","q2 guess","q4 guess","h guess"), 
         col = rep(c("black","red","blue","green"),2), 
         lty = c(1,1,1,1,2,2,2,2), 
         lwd = c(2,2,2,2,1,1,1,1),ncol=2)
  }
  par(oldpar)
}