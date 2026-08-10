module settings_module
  use, intrinsic :: iso_fortran_env, only: dp => real64 
  implicit none

  integer, parameter :: seed = 514, imax = 75, tmax = 33, nmem = 200

  real(dp), parameter :: &
    tau = 1.0_dp, sgm = 0.05_dp, dx = 1.0_dp, &
    d0 = 1.0_dp, d2 = 0.5_dp, d4 = 0.125_dp, &
    k0 = 1.0_dp, k2 = 3.0_dp, k4 = 7.0_dp, kf = 2.0_dp, &
    c0 = 1.0_dp, c2 = -c0 / 3.0_dp, c4 = -c0 / 7.0_dp, cf = 0.5_dp * c0, &
    gma0 = c0 * tau / dx, gma2 = c2 * tau / dx, gma4 = c4 * tau / dx, &
    gmaf = cf * tau / dx, &
    eps0 = 1.0_dp, eps2 = 1.0_dp, eps4 = -1.0_dp

end module settings_module
