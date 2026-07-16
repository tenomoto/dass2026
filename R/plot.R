plot_waves <- function(x, q0, q2, q4, f, i) {
  qr <- range(c(q0, q2, q4, h, f))
  h <- calc_h(q0, q2, q4)
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

