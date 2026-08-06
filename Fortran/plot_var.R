source("../R/settings.R")
source("../R/eqocean.R")
source("../R/plot.R")


read_state <- function(fname) {
  con <- file(fname, "rb")
    q0 <- matrix(readBin(con, "double", imax * tmax), imax, tmax)
    q2 <- matrix(readBin(con, "double", imax * tmax), imax, tmax)
    q4 <- matrix(readBin(con, "double", imax * tmax), imax, tmax)
    f <- matrix(readBin(con, "double", imax * tmax), imax, tmax)
  close(con)
  list(q0 = q0, q2 = q2, q4 = q4, f = f)
}

state <- read_state("state_var.dat")
true_state <- read_state("true_state.dat")

png("var.png", width = 7, height = 7, units = "in", res = 300)
plot_waves(x, state, true_state)
dev.off()

