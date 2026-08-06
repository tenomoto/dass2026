module eqocean_module
  use, intrinsic :: iso_fortran_env, only: dp => real64 
  use settings_module, only: imax
  implicit none

contains

  function eo_forward(q, df, sgm, gma, tau) result(qout)
    real(dp), dimension(:), intent(in) :: q, df
    real(dp), intent(in) :: sgm, gma, tau
    real(dp), dimension(size(q)) :: qout

    qout = (q - 0.5_dp * gma * (cshift(q, 1) - cshift(q, -1)) + &
       0.5_dp * gma ** 2 *(cshift(q, 1) - 2.0_dp * q + cshift(q, -1)) + &
       tau * df) / (1.0_dp + sgm * tau)

  end function eo_forward

  function eo_adjoint(p, dh, sgm, gma, eps, tau) result(pout)
    real(dp), dimension(:), intent(in) :: p, dh
    real(dp), intent(in) :: sgm, gma, eps, tau
    real(dp), dimension(size(p)) :: pout

    pout = (p - 0.5_dp * gma * (cshift(p, -1) - cshift(p, 1)) + &
       0.5_dp * gma ** 2 * (cshift(p, -1) - 2.0_dp * p + cshift(p, 1)) + &
       eps * dh) / (1.0_dp + sgm * tau)
  end function eo_adjoint

  function eo_calc_h(q0, q2, q4) result(h)
    real(dp), dimension(:), intent(in) :: q0, q2, q4
    real(dp), dimension(size(q0)) :: h

    h = q0 + q2 - q4
  
  end function eo_calc_h

end module eqocean_module
