module run_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use settings_module, only: imax, tmax, nmem, &
    tau, sgm, dx, d0, d2, d4, k0, k2, k4, kf, c0, c2, c4, cf, &
    gma0, gma2, gma4, gmaf, eps0, eps2, eps4
  use eqocean_module, only: eo_forward, eo_adjoint, eo_calc_h
  implicit none

  type :: state_type
    real(dp), dimension(imax, tmax) :: q0
    real(dp), dimension(imax, tmax) :: q2
    real(dp), dimension(imax, tmax) :: q4
    real(dp), dimension(imax, tmax) :: f
  end type state_type

  type :: astate_type
    real(dp), dimension(imax) :: p0
    real(dp), dimension(imax) :: p2
    real(dp), dimension(imax) :: p4
    real(dp), dimension(imax, tmax) :: ps
  end type astate_type

  type:: ens_type
    real(dp), dimension(3 * imax) :: xf
    real(dp), dimension(3 * imax, nmem) :: xe
    real(dp), dimension(imax) :: f
    type(state_type) :: mstate
  end type ens_type

contains

function run_forward(q0, q2, q4, f, nt, progf) result(st)
  real(dp), dimension(:), intent(in) :: q0, q2, q4
  real(dp), dimension(:, :), intent(in) :: f
  integer, intent(in) :: nt
  logical, intent(in), optional :: progf
  type(state_type) :: st

  real(dp), dimension(imax) :: df0 = 0.0_dp
  integer :: n

  if (.not. present(progf) .or. progf) then
    st%f(:, 1) = f(:, 1)
    df0(:) = 0.0_dp
    do n = 2, nt
      st%f(:, n) = eo_forward(st%f(:, n - 1), df0, 0.0_dp, gmaf, tau)
    end do
  else
    st%f  = f
  end if
  st%q0(:, 1) = q0(:)
  st%q2(:, 1) = q2(:)
  st%q4(:, 1) = q4(:)
  do n = 2, nt
    st%q0(:, n) = eo_forward(st%q0(:, n-1), d0 * st%f(:, n), sgm, gma0, tau)
    st%q2(:, n) = eo_forward(st%q2(:, n-1), d2 * st%f(:, n), sgm, gma2, tau)
    st%q4(:, n) = eo_forward(st%q4(:, n-1), d4 * st%f(:, n), sgm, gma4, tau)
  end do

end function run_forward

function run_adjoint(dh, ds) result(ast)
  real(dp), dimension(:, :), intent(in) :: dh
  real(dp), dimension(:, :), intent(in), optional :: ds
  type(astate_type) :: ast

  integer :: n

  ast%p0(:) = 0.0_dp ! tmax + 1
  ast%p2(:) = 0.0_dp
  ast%p4(:) = 0.0_dp
  do n = tmax, 1, -1
    ast%p0 = eo_adjoint(ast%p0, dh(:, n), sgm, gma0, eps0, tau)
    ast%p2 = eo_adjoint(ast%p2, dh(:, n), sgm, gma2, eps2, tau)
    ast%p4 = eo_adjoint(ast%p4, dh(:, n), sgm, gma4, eps4, tau)
    if (present(ds) .and. n > 1) then
      ast%ps(:, n) = -tau * (d0 * ast%p0 + d2 * ast%p2 + d4 * ast%p4) - ds(:, n)
    else
      ast%ps(:, n) = 0.0_dp
    end if
  end do

end function run_adjoint

function run_ensemble(xf_in, xe_in, f, nt) result(ens)
  real(dp), dimension(:), intent(in) :: xf_in
  real(dp), dimension(:, :), intent(in) :: xe_in, f
  integer, intent(in) :: nt
  type(ens_type) :: ens

  integer :: mem
  real(dp), dimension(imax) :: q0m, q2m, q4m, q0, q2, q4
  real(dp), dimension(3 * imax) :: xf
  real(dp), dimension(3 * imax, nmem) :: xe
  type(state_type) :: state, mstate

  xf(:) = xf_in(:)
  xe(:, :) = xe_in(:, :)
  q0m = xf(1:imax)
  q2m = xf(imax + 1:2 * imax)
  q4m = xf(2 * imax + 1: 3 * imax)
  do mem = 1, nmem
    q0 = q0m + xe(1:imax, mem)
    q2 = q2m + xe(imax + 1:2 * imax, mem)
    q4 = q4m + xe(2 * imax + 1:3 * imax, mem)
    state = run_forward(q0, q2, q4, f, nt)
    xe(1:imax, mem) = state%q0(:, nt)
    xe(imax + 1:2 * imax, mem) = state%q2(:, nt)
    xe(2 * imax + 1:3*imax, mem) = state%q4(:, nt)
    mstate%q0 = mstate%q0 + state%q0
    mstate%q2 = mstate%q2 + state%q2
    mstate%q4 = mstate%q4 + state%q4
    mstate%f = mstate%f + state%f
  end do
  xf = sum(xe, dim = 2) / nmem
  do mem = 1, nmem
    xe(:, mem) = xe(:, mem) - xf
  end do
  mstate%q0 = mstate%q0 / nmem
  mstate%q2 = mstate%q2 / nmem
  mstate%q4 = mstate%q4 / nmem
  mstate%f = mstate%f / nmem

  ens%xf = xf
  ens%xe = xe
  ens%f = mstate%f(:, tmax)

end function run_ensemble

function run_calc_cost(dh, ds, sh, ss) result(cost)
  real(dp), dimension(:, :), intent(in) :: dh, ds
  real(dp), intent(in) :: sh, ss
  real(dp) :: cost

  cost = 0.5_dp * (sum(dh ** 2) / sh ** 2 + sum(ds ** 2) / ss ** 2)

end function run_calc_cost

subroutine run_print_diag(epoch, cost, gnorm, message) 
  integer, intent(in) :: epoch
  real(dp), intent(in) :: cost, gnorm
  character(len = *), intent(in), optional :: message

  if (present(message)) then
    print *, "epoch: ", epoch, "cost=", cost, "gnorm=", gnorm, " ", message
  else
    print *, "epoch: ", epoch, "cost=", cost, "gnorm=", gnorm
  end if

end subroutine run_print_diag

subroutine run_print_mse(state, tstate)
  type(state_type), intent(in) :: state, tstate

  real(dp), dimension(imax, tmax) :: ha, ht

  ha = eo_calc_h(state%q0, state%q2, state%q4)
  ht = eo_calc_h(tstate%q0, tstate%q2, tstate%q4)
  print *, "MSE = ", sum((ha(:, tmax) - ht(:, tmax)) ** 2)

end subroutine run_print_mse

end module run_module
