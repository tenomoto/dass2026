source("settings.R")
source("plot.R")

if (!exists("hist")) hist <- readRDS("free.rds")

tmax <- ncol(hist$q0)

png("free.png", 800, 1600)
oldpar <- par(mfcol = c(tmax %/% 4, 4), mar = c(2, 2, 1, 1))

for (t in 1:tmax) {
  plot_waves(x, hist$q0[, t], hist$q2[, t], hist$q4[, t], hist$f[, t], t)
}
dev.off()

par(oldpar)
