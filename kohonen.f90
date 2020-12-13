module kohonen_module
  implicit none

  integer, parameter :: kohonen_eckert = 0, kohonen_gaussian = 1

  integer, private :: imax, jmax, kmax, tmax, pmax
  real, private :: sigma0, tau0, ptau, smin
  integer, private :: p = 1
  real, dimension(:, :, :, :), allocatable, private :: delta
  integer, dimension(2), private :: ns

  real, dimension(:, :), allocatable, public :: kohonen_drms
  real, dimension(:, :, :), allocatable, public :: kohonen_weight
  logical, public :: kohonen_debug = .false.

  public :: kohonen_Init, kohonen_Normalize, kohonen_Learn, kohonen_Elect_unit


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


  subroutine kohonen_Init(kx, ix, jx, tx, px, x, y, s0, t0, w0, tau, smn, lhex)
    integer, intent(in) :: kx, ix, jx, tx, px
    real, intent(out), dimension(:,:), allocatable :: x
    real, intent(out), dimension(:), allocatable :: y
    real, intent(in), optional :: s0, t0, w0, tau, smn
    logical, intent(in), optional :: lhex

    real, parameter :: sfact = 0.2, tau00 = 0.25
    real :: w00 = 1.0e-5
    real :: wmean, wsdev
    real, dimension(:), allocatable :: w
    integer :: i, j, ii, jj
    logical :: lhexagonal = .true.

    ! Initialize parameters and allocate arrays
    kmax = kx; imax = ix; jmax = jx
    tmax = tx; pmax = px
    if ( kohonen_debug ) then
      print *, "kmax =", kmax, "imax =", imax, " jmax = ", jmax, &
        "tmax =", tmax, " pmax = ", pmax
    end if
    allocate( delta(imax, jmax, imax, jmax), &
      kohonen_weight(kmax, imax, jmax), &
      kohonen_drms(imax, jmax))
    allocate( x(imax, jmax), y(jmax), w(imax * jmax * kmax) )

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
    kohonen_weight(:, :, :) = reshape( w, (/kmax, imax, jmax/) )

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
      sigma0 = 0.5 * sqrt( maxval(x) ** 2 + maxval(y) ** 2)
    end if
    if (present(t0)) then
      tau0 = t0
    else
      tau0 = tau00
    end if
    if ( present(tau) ) then
      ptau = tau
    else
      ptau = real(pmax)
    end if
    if ( present(smn) ) then
      smin = smn
    else
      smin = sfact * sigma0
    end if
    if ( kohonen_debug ) then
      print *, "sigma0 =", sigma0, " tau0 =", tau0, " ptau =", ptau, " smin=", smin
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


  subroutine kohonen_Elect_unit( z, u, d_p )
    real, dimension(:), intent(in) :: z
    integer, dimension(2), intent(out) :: u 
    real, intent(out) :: d_p

    integer :: i, j

    do j = 1, jmax
      do i = 1, imax
        kohonen_drms(i, j) = sqrt( sum( (kohonen_weight(:, i, j) - z(:)) ** 2 ) )
      end do
    end do
    d_p = minval(kohonen_drms)
    u = minloc(kohonen_drms)
    if ( kohonen_debug )  then
      print *, "Learning step p =", p, "/", pmax
      print *, "elected unit =", u, " d_p = ", d_p
    end if

  end subroutine kohonen_Elect_unit


  subroutine Adjust_weight_Ekert( z )
    real, dimension(:), intent(in) :: z

    integer :: i, j
    real :: denom, sigma, tau, wmean, wsdev
    real, dimension(:), allocatable :: w

    denom = 1.0 / (1.0 + p / ptau)
    sigma = sigma0 * denom
    tau = tau0 * denom
    do j = 1, jmax
      do i = 1, imax
        if ( delta(i, j, ns(1), ns(2)) <= sigma ) then
          kohonen_weight(:, i, j) = kohonen_weight(:, i, j) + &
            tau * (z(:) - kohonen_weight(:, i, j))
        end if
      end do
    end do
    allocate( w(imax * jmax * kmax) )
    w = pack( kohonen_weight, .true. )
    call kohonen_Normalize( w,  wmean, wsdev )
    kohonen_weight(:, :, :) = reshape( w, (/kmax, imax, jmax/) )
    if ( kohonen_debug )  then
      print *, "sigma =", sigma, " tau =", tau
    end if

  end subroutine Adjust_weight_Ekert


  subroutine Adjust_weight_Gaussian( z )
    real, dimension(:), intent(in) :: z

    integer :: i, j
    real :: sigma
    real, dimension(:, :), allocatable :: h

    allocate(h(imax, jmax))
    sigma = max( sigma0 * exp( -p / ptau), smin )
    h(:, :) = exp( -0.5 * (delta(:, :, ns(1), ns(1)) / sigma) ** 2 )
    do j = 1, jmax
      do i = 1, imax
        kohonen_weight(:, i, j) = kohonen_weight(:, i, j) + &
            h(i, j) * (z(:) - kohonen_weight(:, i, j))
      end do
    end do
    if ( kohonen_debug )  then
      print *, "sigma =", sigma, " hmax=", maxval(h), " hmin=", minval(h)
    end if
    deallocate(h)

  end subroutine Adjust_weight_Gaussian


  subroutine kohonen_Learn( z, d_p, m )
    real, dimension(:), intent(in) :: z
    real, intent(out) :: d_p
    integer :: m

    call kohonen_Elect_unit( z, ns, d_p )
    if ( m == kohonen_eckert ) then
      call Adjust_weight_Ekert( z )
    else
      call Adjust_weight_Gaussian( z )
    end if

    p = p + 1

  end subroutine kohonen_Learn


end module kohonen_module
