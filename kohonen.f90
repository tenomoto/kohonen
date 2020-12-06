module kohonen_module
  implicit none

!  integer, parameter
!    imax = 8, jmax = 8, kmax = 10, lmax = 10, tmax = 2192, pmax = 5.0 * tmax
!  real, parameter :: sigma0 = 5.66, tau0 = 0.25
  real, parameter :: tau00 = 0.25

  integer, private :: imax, jmax, kmax, lmax, tmax, pmax
  real, private :: sigma0, tau0
  integer, private :: p = 1
  real, dimension(:, :, :, :), allocatable, private :: delta
  integer, dimension(2), private :: ns

  real, dimension(:, :), allocatable, public :: kohonen_drms
  real, dimension(:, :, :, :), allocatable, public :: kohonen_weight
  logical, public :: kohonen_debug = .false.

  public :: kohonen_Init, kohonen_Normalize, kohonen_Learn


contains
  

  subroutine kohonen_Normalize(x, xmean, xsdev)
    real, dimension(:), intent(inout) :: x 
    real, intent(out) :: xmean, xsdev

    integer :: nmax

    nmax = size( x(:) )
    xmean = sum( x(:) ) / nmax
    xsdev = sqrt( (sum( x(:) ** 2 )  - nmax * xmean **2) / (nmax - 1.0) )
    x = (x(:) - xmean) / xsdev

  end subroutine kohonen_Normalize


  subroutine kohonen_Init(kx, lx, ix, jx, tx, px, x, y, s0, t0, w0, lhex)
    integer, intent(in) :: kx, lx, ix, jx, tx, px
    real, intent(in), optional :: s0, t0, w0
    real, intent(out), dimension(:,:), allocatable :: x
    real, intent(out), dimension(:), allocatable :: y
    logical, intent(in), optional :: lhex

    real :: w00 = 1.0e-5
    real :: wmean, wsdev
    real, dimension(:), allocatable :: w
    integer :: i, j, ii, jj
    logical :: lhexagonal = .true.

    ! Initialize parameters and allocate arrays
    kmax = kx; lmax = lx
    imax = ix; jmax = jx
    tmax = tx; pmax = px
    if ( kohonen_debug ) then
      print *, "kmax =", kmax, "lmax =", lmax, &
        "imax =", imax, " jmax = ", jmax, &
        "tmax =", tmax, " pmax = ", pmax
    end if
    allocate( delta(imax, jmax, imax, jmax), &
      kohonen_weight(kmax, lmax, imax, jmax), &
      kohonen_drms(imax, jmax))
    allocate( x(imax, jmax), y(jmax), w(imax * jmax * kmax * lmax) )

    ! Initialize kohonen_weight with small random numbers
    call random_number( w )
    call kohonen_Normalize( w,  wmean, wsdev )
    if ( present(w0) ) then
      w00 = w0
    end if
    w(:) = w00 * w(:)
    if ( kohonen_debug ) then
      print *, "weight max =", maxval(w), " min = ", minval(w), &
        " amplitude =", w00
    end if
    kohonen_weight(:, :, :, :) = reshape( w, (/kmax, lmax, imax, jmax/) )

    if (present(lhex)) then
      lhexagonal = lhex
    end if
    ! rectangular unit layout
    do j = 1, jmax 
      do i = 1, imax 
        x(i, j) = real(i) - 1.0
      end do
      y(j) = real(j) - 1.0
    end do
    if ( lhexagonal ) then
      ! shift right for odd columns
      x(:, 1::2) = x(:, 1::2) + 0.5
      ! height of equilateral triangle
      y(:) = y(:) * 0.5 * sqrt(3.0)
    end if
    if ( kohonen_debug ) then
      do j = 1, jmax
        print *, "x =", x(:, j), " y = ", y(j)
      end do
    end if
    if ( kohonen_debug ) then
      if ( lhexagonal ) then
        print *, "hexagonal layout"
      else
        print *, "rectangular layout"
      end if
    end if

    if (present(s0)) then
      sigma0 = s0
    else
      ! half the diagonal of the domain
      print *, maxval(x), maxval(y)
      sigma0 = 0.5 * sqrt( maxval(x) ** 2 + maxval(y) ** 2)
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
            delta(i, j, ii, jj) = sqrt( (x(i, j) - x(ii, jj)) ** 2 + (y(j) - y(jj)) ** 2 )
          end do
        end do
      end do
    end do

    deallocate(w)

  end subroutine kohonen_Init


  subroutine Elect_unit( z, d_p )
    real, dimension(:, :), intent(in) :: z
    real, intent(out) :: d_p

    integer :: i, j

    ! Elect unit
    do j = 1, jmax
      do i = 1, imax
        kohonen_drms(i, j) = sqrt( sum( (kohonen_weight(:, :, i, j) - z(:, :))**2 ) )
      end do
    end do
    d_p = minval(kohonen_drms)
    ns = minloc(kohonen_drms)
    if ( kohonen_debug )  then
      print *, "Learning step p =", p, "/", pmax
      print *, "elected unit =", ns, " d_p = ", d_p
    end if

  end subroutine Elect_unit

  subroutine Adjust_weight_Ekert( z )
    real, dimension(:, :), intent(in) :: z

    integer :: i, j
    real :: denom, sigma, tau, wmean, wsdev
    real, dimension(:), allocatable :: w

    denom = 1.0 / (1.0 + p / pmax)
    sigma = sigma0 * denom
    tau = tau0 * denom
    do j = 1, jmax
      do i = 1, imax
        if ( delta(i, j, ns(1), ns(2)) <= sigma ) then
          kohonen_weight(:, :, i, j) = kohonen_weight(:, :, i, j) + &
            tau * (z(:, :) - kohonen_weight(:, :, i, j))
        end if
      end do
    end do
    allocate( w(imax * jmax * kmax * lmax) )
    w = pack( kohonen_weight, .true. )
    call kohonen_Normalize( w,  wmean, wsdev )
    kohonen_weight(:, :, :, :) = reshape( w, (/kmax, lmax, imax, jmax/) )
    if ( kohonen_debug )  then
      print *, "sigma =", sigma, " tau =", tau
    end if

  end subroutine Adjust_weight_Ekert


  subroutine kohonen_Learn( z, d_p )
    real, dimension(:, :), intent(in) :: z
    real, intent(out) :: d_p

    call Elect_unit( z, d_p )

    call Adjust_weight_Ekert( z )

    p = p + 1

  end subroutine kohonen_Learn


end module kohonen_module
