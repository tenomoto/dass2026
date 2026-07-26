module random_module
  use, intrinsic :: iso_fortran_env, only: dp => real64

contains

  subroutine random_set_seed(seed)
    integer :: seed
    integer :: n
    integer, allocatable :: s(:)

    call random_seed(size=n)
    print *, "seed size=", n
    allocate(s(n))
    s = seed
    call random_seed(put=s)
    call random_seed(get=s)
    print *, "seed=", s
    deallocate(s)

  end subroutine random_set_seed

  function random_uniform(n) result(u)
    integer, intent(in) :: n
    real(dp), allocatable :: u(:)

    allocate(u(n))
    call random_number(u)
    u = 1 - u

  end function random_uniform

  function random_normal(n, mean, sd) result(x)
    integer, intent(in) :: n 
    real(dp), intent(in), optional :: mean, sd
    real(dp), allocatable :: x(:)

    real(dp), allocatable :: u1(:), u2(:)
    real(dp) :: tau

    tau = 2.0_dp * acos(-1.0_dp)
    allocate(u1(n), u2(n), x(n))
    u1 = random_uniform(n)
    u2 = random_uniform(n)
    x = sqrt(-2.0_dp * log(u1)) * cos(tau * u2)
    if (present(sd)) x = sd * x
    if (present(mean)) x = mean + x
    deallocate(u1, u2)

  end function random_normal

  subroutine random_test()
    integer, parameter :: n = 1000000, seed = 514
    real(dp), parameter :: mu = 1.0_dp, sig = 4.0_dp
    real(dp), allocatable :: u(:), x(:)
    real(dp) :: xmean, xstd

    call random_set_seed(seed)
    u = random_uniform(n)
    print *, "uniform: nzero=", count(u==0) 
    x = random_normal(n)
    xmean = sum(x) / size(x)
    xstd = sqrt(sum((x - xmean) ** 2) / size(x))
    print *, "normal: xmean=", xmean, " xstd=", xstd
    x = random_normal(n, mu, sig)
    xmean = sum(x) / size(x)
    xstd = sqrt(sum((x - xmean) ** 2) / size(x))
    print *, "normal mean=", mu, " sig=", sig, ": xmean=", xmean, " xstd=", xstd 
    open(unit=51, file="rand.txt", status="replace")
    write(unit=51, fmt=*) x
    close(unit=51)

  end subroutine random_test

end module random_module
  
