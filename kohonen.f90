module kohonen_module
  implicit none

!  integer, parameter
!    imax = 8, jmax = 8, kmax = 10, lmax = 10, tmax = 2192, pmax = 5.0 * tmax
!  real, parameter :: sigma0 = 5.66, tau0 = 0.25
  real, parameter :: tau00 = 0.25

  integer, private :: imax, jmax, kmax, lmax, tmax, pmax
  real, private :: sigma0, tau0
  integer, private :: p = 1
  real, dimension(:), allocatable, private :: x, y
  real, dimension(:, :, :, :), allocatable, private :: delta
  integer, dimension(2), private :: ns

  real, dimension(:, :), allocatable, public :: kohonen_drms
  real, dimension(:, :, :, :), allocatable, public :: kohonen_weight
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
    integer :: nmax

    nmax = size( kohonen_weight )
    wmean = sum( kohonen_weight(:, :, :, :) ) / nmax
    wsdev = sqrt( (sum( kohonen_weight(:, :, :, :) ** 2 ) &
      - nmax * wmean **2) / (nmax - 1.0) )
    kohonen_weight = (kohonen_weight(:, :, :, :) - wmean) / wsdev

  end subroutine Normalize_weight


  subroutine kohonen_Init(kx, lx, ix, jx, tx, px, s0, t0, lhex)
    integer, intent(in) :: kx, lx, ix, jx, tx, px
    real, intent(in), optional :: s0, t0
    logical, intent(in), optional :: lhex

    integer :: i, j, ii, jj

    ! Initialize parameters and allocate arrays
    kmax = kx; lmax = lx
    imax = ix; jmax = jx
    tmax = tx; pmax = px
    if ( kohonen_debug ) then
      print *, "kmax =", kmax, "lmax =", lmax, &
        "imax =", imax, " jmax = ", jmax, &
        "tmax =", tmax, " pmax = ", pmax
    end if
    allocate( x(imax), y(jmax), &
      delta(imax, jmax, imax, jmax), &
      kohonen_weight(kmax, lmax, imax, jmax), &
      kohonen_drms(imax, jmax))

    ! Initialize kohonen_weight with small random numbers
    call random_number(kohonen_weight)
    call Normalize_weight()

    ! rectangular unit layout
    do i = 1, imax 
      x(i) = real(i) - 1.0
    end do
    do j = 1, jmax 
      y(j) = real(j) - 1.0
    end do
    ! hexagonal unit layout
    if (present(lhex)) then
      if ( lhex ) then
        do i = 1, imax, 2 
          x(i) = x(i) + 0.5
        end do
        do j = 1, jmax 
          y(j) = y(j) * 0.5 * sqrt(3.0)
        end do
      end if
      if ( kohonen_debug ) then
        print *, "hexagonal layout"
      end if
    else
      print *, "rectangular layout"
    end if
    if ( kohonen_debug ) then
      print *, "x =", x
      print *, "y =", y
    end if

    if (present(s0)) then
      sigma0 = s0
    else
      ! half the diagonal of the domain
      sigma0 = 0.5 * sqrt( (maxval(x) - minval(x))**2 + (maxval(y) - minval(y)) ** 2)
    end if
    if (present(t0)) then
      tau0 = t0
    else
      tau0 = tau00
    end if
    if ( kohonen_debug ) then
      print *, "sigma0 =", sigma0, " tau0 =", tau0
    end if

    ! Calculate delta as Cartesian distance between units
    do jj = 1, jmax
      do ii = 1, imax
        do j = 1, jmax
          do i = 1, imax
            delta(i, j, ii, jj) = sqrt( (x(i) - x(ii)) ** 2 + (y(j) - y(jj)) ** 2 )
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
        kohonen_drms(i, j) = sqrt( sum( (kohonen_weight(:, :, i, j) - z(:, :))**2 ) )
      end do
    end do
    dp = minval(kohonen_drms)
    ns = minloc(kohonen_drms)

    ! Adjust kohonen_weight 
    call Calc_sigma_tau( p, sigma, tau )
    do j = 1, jmax
      do i = 1, imax
        if ( delta(i, j, ns(1), ns(2)) <= sigma ) then
          kohonen_weight(i, j, :, :) = kohonen_weight(i, j, :, :) + &
            tau * (z(:, :) - kohonen_weight(:, :, i, j))
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
