module kohonen_module
  implicit none

!  integer, parameter
!    imax = 8, jmax = 8, kmax = 10, lmax = 10, tmax = 2192, pmax = 5.0 * tmax
!  real, parameter :: sigma0 = 5.66, tau0 = 0.25

  integer, private :: imax, jmax, kmax, lmax
  real, private :: tmax, pmax, sigma0, tau0, nmax
  integer, private :: p = 1
  real, dimension(:, :), allocatable, private :: drms
  real, dimension(:, :, :, :), allocatable, private :: delta, weight
  integer, dimension(2), private :: ns

  logical, public :: kohonen_debug = .false.

  public :: kohonen_Init, kohonen_Learn

contains
  
  subroutine Calc_sigma_tau( p, sigma, tau )
    integer, intent(in) :: p
    real, intent(out) :: sigma, tau

    real :: denom

    denom = 1.0 / (1.0 + p / pmax)
    sigma = sigma0 * denom
    tau = tau0 * denom

  end subroutine Calc_sigma_tau


  subroutine Normalize_weight()
    real :: wmean, wsdev

    wmean = sum( weight(:, :, :, :) ) / nmax
    wsdev = sqrt( (sum( weight(:, :, :, :)**2 ) - nmax * wmean **2) / (nmax - 1.0) )
    weight = (weight(:, :, :, :) - wmean) / wsdev

  end subroutine Normalize_weight


  subroutine kohonen_Init(kx, lx, ix, jx, tx, px, s0, t0)

    integer, intent(in) :: kx, lx, ix, jx, tx, px
    real, intent(in) :: s0, t0

    integer :: i, j, ii, jj

    ! Initialize parameters and allocate arrays
    kmax = kx; lmax = lx
    imax = ix; jmax = jx
    tmax = tx; pmax = px
    nmax = imax * jmax * kmax * lmax
    sigma0 = s0; tau0 = t0
    allocate( &
      delta(imax, jmax, imax, jmax), &
      weight(kmax, lmax, imax, jmax), &
      drms(imax, jmax))

    ! Initialize weight with small random numbers
    call random_number(weight)
    call Normalize_weight()

    ! Calculate delta as Cartesian distance between units
    do jj = 1, jmax
      do ii = 1, imax
        do j = 1, jmax
          do i = 1, imax
            delta(i, j, ii, jj) = sqrt(real( (i - ii) ** 2 + (j - jj) ** 2 ))
          end do
        end do
      end do
    end do

  end subroutine kohonen_Init


  subroutine kohonen_Learn( z, dp )

    real, dimension(:, :), intent(in) :: z
    real, intent(out) :: dp

    integer :: i, j
    real :: sigma, tau

    ! Elect unit
    do j = 1, jmax
      do i = 1, imax
        drms(i, j) = sqrt( sum( (weight(:, :, i, j) - z(:, :))**2 ) )
      end do
    end do
    dp = minval(drms)
    ns = minloc(drms)

    ! Adjust weight 
    call Calc_sigma_tau( p, sigma, tau )
    do j = 1, jmax
      do i = 1, imax
        if ( delta(i, j, ns(1), ns(2)) <= sigma ) then
          weight(i, j, :, :) = weight(i, j, :, :) + &
            tau * (z(:, :) - weight(:, :, i, j))
        end if
      end do
    end do
    call Normalize_weight()

    if ( kohonen_debug )  then
      print *, "Learning step p =", p, "/", pmax
      print *, "elected unit =", ns
    end if

    p = p + 1

  end subroutine kohonen_Learn


end module kohonen_module
