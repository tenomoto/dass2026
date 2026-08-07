program enkf
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use matrix_module, only: matrix_inv
  use random_module, only: random_set_seed, random_normal
  use settings_module, only: k0, k2, k4, kf, imax, tmax, nmem, seed
  use eqocean_module, only: eo_calc_h
  use run_module, only: state_type, ens_type, &
    run_forward, run_ensemble, run_calc_cost, &
    run_print_diag, run_print_mse, run_save_state
  implicit none

  integer, parameter :: ntobs = 3
  integer, parameter, dimension(ntobs) :: tobs = [0, 16, 32] + 1
  real(dp), parameter :: pi = acos(-1.0_dp), &
    s_ens = 0.01_dp, s_obs = 0.0_dp

  integer :: i, mem, n, k, n_next, nt
  integer, dimension(ntobs - 1) :: dtobs
  real(dp), dimension(imax) :: x, q0, q2, q4
  real(dp), dimension(3 * imax) :: xf, xf0
  real(dp), dimension(imax, ntobs) :: hobs
  real(dp), dimension(imax, tmax) :: h, htrue,  &
    f = 0.0_dp, dh = 0.0_dp, ds = 0.0_dp, &
    q0_hist, q2_hist, q4_hist, f_hist
  real(dp), dimension(imax, imax) :: rmat
  real(dp), dimension(3 * imax, nmem) :: xe
  type(state_type) :: true_state, state
  type(ens_type) :: estate, efcst

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
  hobs = htrue(:, tobs) + reshape(random_normal(imax * ntobs, 0.0_dp, s_obs), [imax, ntobs])

  q0 = sin(k0 * x)
  q2 = sin(k2 * x)
  q4 = sin(k4 * x)
  f = 0.0_dp
  f(:, 1) = 0.1 * cos(kf * x)

  rmat = 0.0_dp
  do i = 1, imax
    rmat(i, i) = s_obs**2
  end do
  xf(1:imax) = q0
  xf(imax + 1:2 * imax) = q2
  xf(2 * imax + 1:3 * imax) = q4
  xf0 = xf
  xe = reshape(random_normal(3 * imax * nmem, 0.0_dp, s_ens), [3 * imax, nmem])
  f = true_state%f

  dtobs = tobs(2:ntobs) - tobs(1:ntobs - 1)
  n = 1
  do k = 1, ntobs
    print *, "k=", k,"n=", n, tobs(k)
    if (n == tobs(k)) then
      print *, "update ensemble"
      estate = ensemble_update(xf, xe, hobs(:, k), rmat)
      xf = estate%xf
      xe = estate%xe
      q0_hist(:, n) = xf(1:imax)
      q2_hist(:, n) = xf(imax + 1:2 * imax)
      q4_hist(:, n) = xf(2 * imax + 1:3*imax)
    end if
    if (k < ntobs) then 
      print *, "run ensemble forecast"
      n_next = n + dtobs(k)
      nt = dtobs(k) + 1
      f =  true_state%f(:, n:n_next)
      efcst = run_ensemble(xf, xe, f, nt)
      xf = efcst%xf
      xe = efcst%xe
      q0_hist(:, (n + 1):n_next) = efcst%mstate%q0(:, 2:nt)
      q2_hist(:, (n + 1):n_next) = efcst%mstate%q2(:, 2:nt)
      q4_hist(:, (n + 1):n_next) = efcst%mstate%q4(:, 2:nt)
      f_hist(:, (n + 1):n_next) = efcst%mstate%f(:, 2:nt)
      n = n_next
    end if
  end do
  state%q0 = q0_hist
  state%q2 = q2_hist
  state%q4 = q4_hist
  state%f = f_hist
  call run_print_mse(state, true_state)
  call run_save_state("state_enkf.dat", state)

contains

  function ensemble_update(xf_in, xe_in, yo, rmat, infl_in, std_in) result(ens_state)
    real(dp), dimension(:), intent(in) :: xf_in
    real(dp), dimension(:, :), intent(in) :: xe_in
    real(dp), dimension(:), intent(in) :: yo
    real(dp), dimension(:, :), intent(in) :: rmat
    real(dp), optional :: infl_in, std_in
    type(ens_type) :: ens_state

    integer :: omax, mem
    real(dp), dimension(size(yo)) :: hxm
    real(dp), dimension(size(yo), nmem) :: hxe, ye
    real(dp), dimension(3 * imax, size(yo)) :: kmat
    real(dp), dimension(imax) :: x, q0, q2, q4
    real(dp), dimension(3 * imax) :: xf
    real(dp), dimension(3 * imax, nmem) :: xe
    real(dp) :: infl, std

    if (.not. present(infl_in)) then
      infl = 0.1_dp
    else
      infl = infl_in
    end if
    if (.not. present(std_in)) then
      std = 0.01_dp
    else
      std = std_in
    end if

    xf = xf_in
    xe = xe_in

    omax = size(yo)
    hxm = 0.0_dp
    xe = sqrt(1.0_dp + infl) * xe
    do mem = 1, nmem
      q0 = xf(1:imax) + xe(1:imax, mem)
      q2 = xf(imax + 1:2 * imax) + xe(imax + 1:2 * imax, mem)
      q4 = xf(2 * imax + 1:3 * imax) + xe(2 * imax + 1:3 * imax, mem)
      hxe(:, mem) = eo_calc_h(q0, q2, q4)
      hxm = hxm + hxe(:, mem)
    end do
    hxm = hxm / nmem
    do mem = 1, nmem
      hxe(:, mem) = hxe(:, mem) - hxm
    end do
    ! K = X^f*Y^T*(Y*Y^T + (N-1)*R)^{-1}
    kmat = matmul(xe, matmul(transpose(hxe),  &
      matrix_inv(matmul(hxe, transpose(hxe)) + (nmem - 1) * rmat)))
    xf = xf + matmul(kmat, yo - hxm)
    ye = reshape(random_normal(omax * nmem, 0.0_dp, std), [omax, nmem])
    xe = xe + matmul(kmat, ye - hxe)
    ens_state%xf = xf
    ens_state%xe = xe

  end function ensemble_update

end program enkf
