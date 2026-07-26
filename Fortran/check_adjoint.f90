program check_adjoint
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: random_set_seed, random_normal
  use settings_module, only: seed, imax, sgm, gma0, tau
  use eqocean_module, only: eo_forward, eo_adjoint
  implicit none

  integer, parameter :: tmax = 33
  real(dp), parameter :: eps = 0.0_dp

  real(dp), dimension(imax) :: qi, df, dh, eta, p
  real(dp), dimension(imax, tmax) :: q
  real(dp) :: lhs, rhs
  integer :: i

  call random_set_seed(seed)

  qi(:) = random_normal(imax)
  df(:) = 0.0_dp ! tlm has not forcing
  dh(:) = 0.0_dp ! no obs
  eta(:) = 0.0_dp

  q(:, :) = 0.0_dp
  q(:, 1) = qi(:)
  do i = 2, tmax
    q(:, i) = eo_forward(q(:, i - 1), df, sgm, gma0, tau)
  end do

  lhs = sum(q(:, tmax) * q(:, tmax))
  print *, "LHS:", lhs
  p = q(:, tmax)
  do i = tmax - 1, 1, -1
    p = eo_adjoint(p, dh, sgm, gma0, eps, tau)
  end do
  rhs = sum(qi * p)
  print *, "RHS:", rhs
  print *, "LHS - RHS:", lhs - rhs

end program check_adjoint

 
