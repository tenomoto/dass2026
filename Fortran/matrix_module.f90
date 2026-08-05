module matrix_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

contains

! https://fortranwiki.org/fortran/show/Matrix+inversion
  function matrix_inv(amat) result(ainv)
    real(dp), dimension(:, :), intent(in) :: amat
    real(dp), dimension(size(amat, 1), size(amat, 2)) :: ainv

    real(dp), dimension(size(amat, 1)) :: work
    integer, dimension(size(amat, 1)) :: ipiv
    integer :: n, info

    external dgetrf, dgetri

    ainv = amat
    n = size(amat, 1)

    call dgetrf(n, n, ainv, n, ipiv, info)
    if (info /= 0) then
      stop "Error in inv: amat is singular."
    end if

    call dgetri(n, ainv, n, ipiv, work,n, info)
    if (info /= 0) then
      stop "Error in inv: matrix inversion failed."
    end if

  end function matrix_inv

end module matrix_module

