program var
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: random_set_seed, random_normal
  use settings_module, only: k0, k2, k4, kf, imax, tmax, seed
  use eqocean_module, only: eo_calc_h
  use run_module, only: state_type, astate_type, &
    run_forward, run_adjoint, run_calc_cost, &
    run_print_diag, run_print_mse, run_save_state
  implicit none

  integer, parameter :: ntobs = 3, niter = 1000, fcg = 1
  integer, parameter, dimension(ntobs) :: tobs = [0, 16, 32] + 1
  real(dp), parameter :: pi = acos(-1.0_dp), &
    sobs = 0.0_dp, smod = 0.01_dp, &
    ftol = 1.0e-15_dp, ctol = 1.0e-15_dp, gtol = 1.0e-15_dp, &
    lr = 0.1_dp
  real(dp), parameter :: sh = 1.0_dp, ss = 1.0_dp

  integer :: i, ep, n
  real(dp) :: ognorm = 1.0e8_dp, gnorm, bta, ocost = 1.0e8_dp, cost, dcost
  real(dp), dimension(imax) :: x, q0, q2, q4
  real(dp), dimension(3 * imax) :: d = 0.0_dp, xf, g, xa
  real(dp), dimension(imax, tmax) :: h, htrue, f = 0.0_dp, dh = 0.0_dp, ds = 0.0_dp
  real(dp), dimension(imax, ntobs) :: hobs
  type(state_type) :: true_state, state
  type(astate_type) :: adj

  call random_set_seed(seed)

  x = 2 * pi / imax * [(i, i = 0, imax - 1)]
  q0 = cos(k0 * x)
  q2 = cos(k2 * x)
  q4 = cos(k4 * x)
  f = 0.0_dp 
  f(:, 1) = 0.1 * cos(kf * x)

  true_state = run_forward(q0, q2, q4, f, tmax)

  do n = 1, tmax
   htrue(:, n) = eo_calc_h(true_state%q0(:, n), true_state%q2(:, n), true_state%q4(:, n))
  end do
  hobs = htrue(:, tobs) + reshape(random_normal(imax * ntobs, 0.0_dp, sobs), [imax, ntobs])

  q0 = sin(k0 * x)
  q2 = sin(k2 * x)
  q4 = sin(k4 * x)
  f = 0.0_dp
  f(:, 1) = 0.1 * cos(kf * x)

  do ep = 1, niter
    state = run_forward(q0, q2, q4, f, tmax)
    do n = 1, tmax
     h(:, n) = eo_calc_h(state%q0(:, n), state%q2(:, n), state%q4(:, n))
    end do
    dh(:, tobs) = hobs - h(:, tobs)
    adj = run_adjoint(dh)
    cost = run_calc_cost(dh, ds, sh, ss)
    xf(1:imax) = q0
    xf(imax + 1:2 * imax) = q2
    xf(2 * imax + 1:3 * imax) = q4
    g(1:imax) = -adj%p0
    g(imax + 1:2 * imax) = -adj%p2
    g(2 * imax + 1:3 * imax) = -adj%p4
    gnorm = sum(g ** 2)
    bta = gnorm / ognorm
    d = -g + fcg * bta * d
    xa = xf + lr * d
    q0 = xa(1:imax)
    q2 = xa(imax + 1:2 * imax)
    q4 = xa(2 * imax + 1:3 * imax)
    dcost = abs(cost - ocost)
    if (mod(ep - 1, 100) == 0) then
      call run_print_diag(ep, cost, gnorm)
    end if
    if (cost < ftol) then
      call run_print_diag(ep, cost, gnorm, "stopping due to ftol")
      exit
    end if
    if (dcost < ctol) then
      call run_print_diag(ep, cost, gnorm, "stopping due to ctol")
      exit
    end if
    if (gnorm < gtol) then
      call run_print_diag(ep, cost, gnorm, "stopping due to gtol")
      exit
    end if
    ocost = cost
    ognorm = gnorm
  end do
  state = run_forward(q0, q2, q4, f, tmax)
  call run_print_mse(state, true_state)

  call run_save_state("state_var.dat", state)
  call run_save_state("true_state.dat", true_state)

end program var
