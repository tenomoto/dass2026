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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage:: Rscript plot.R var|enkf")
}
method <- args[1]
state <- read_state(paste0("state_", method, ".dat"))
true_state <- read_state("true_state.dat")

png(paste0(method, ".png"), width = 7, height = 7, units = "in", res = 300)
plot_waves(x, state, true_state)
dev.off()

