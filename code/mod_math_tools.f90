
!=======================================================================
!  This code was adapted from MathToolsMod.F90                         !
!  in the link below.                                                  !
!  https://github.com/gbonan/CLM-ml_v0                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Bonan G. B., Williams M., Fisher R. A., Oleson K. W. (2014),        !
!  Modeling stomatal conductance in the earth system: linking leaf     !
!  water-use efficiency and water transport along the                  !
!  soil–plant–atmosphere continuum, Geoscientific Model Development,   !
!  7, 2019-2222.
!                                                                      !
!=======================================================================

module mod_math_tools

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Math tools
  !
  ! !USES:
  use mod_param
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: zbrent             ! Use Brent's method to find the root of a function (water use efficiency optimization)
  public :: quadratic          ! Solve a quadratic equation for its two roots
  public :: tridiag            ! Solve a tridiagonal system of equations

  interface
    function xfunc (o2air, co2air, iota, i, x, &
                    universal, atmos, layer, flux) result(f)
    use mod_param
    real(8), intent(in) :: o2air, co2air, iota
    integer, intent(in) :: i
    real(8), intent(in) :: x
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    real(8) :: f
    end function xfunc
  end interface

contains

  !-----------------------------------------------------------------------
  function zbrent (msg, func, o2air, co2air, iota, i, &
                   universal, atmos, layer, flux, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Use Brent's method to find the root of a function, which is known to exist
    ! between xa and xb. The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*)    :: msg        ! String to be printed
    procedure (xfunc)   :: func       ! Function to solve
    real(8), intent(in) :: o2air, co2air
    real(8), intent(in) :: iota
    integer, intent(in) :: i
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    real(8), intent(in) :: xa, xb    ! Minimum and maximum of the variable domain to search
    real(8), intent(in) :: tol       ! Error tolerance
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: itmax = 50          ! Maximum number of iterations
    real(8), parameter :: eps = 1.e-08        ! Relative error tolerance
    integer     :: iter                       ! Iteration loop index
    real(8)     :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    real(8)     :: root
    !---------------------------------------------------------------------

    a = xa
    b = xb
    fa = func (o2air, co2air, iota, i, a, universal, atmos, layer, flux)
    fb = func (o2air, co2air, iota, i, b, universal, atmos, layer, flux)

    if ((fa > 0.0d0 .and. fb > 0.0d0) .or. (fa < 0.0d0 .and. fb < 0.0d0)) then
       write (*,*) 'ERROR!!! zbrent: Root must be bracketed'
       write (*,*) 'called from: ',msg
       write (*,*) xa, fa
       write (*,*) xb, fb
    end if
    c = b
    fc = fb
    iter = 0
    do
!write(*,*) 'iter, b, fb', iter, b, fb
       if (iter == itmax) exit
       iter = iter + 1
       if ((fb > 0.0d0 .and. fc > 0.0d0) .or. (fb < 0.0d0 .and. fc < 0.0d0)) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
       xm = 0.5d0 * (c - b)
       if (abs(xm) <= tol1 .or. fb == 0.0d0) exit
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s = fb / fa
          if (a == c) then
             p = 2.0d0 * xm * s
             q = 1.0d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2.0d0 * xm * q * (q - r) - (b-a) * (r - 1.0d0))
             q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
          end if
          if (p > 0.0d0) q = -q
          p = abs(p)
          if (2.0d0*p < min(3.0d0*xm*q-abs(tol1*q),abs(e*q))) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1,xm)
       end if
       fb = func (o2air, co2air, iota, i, b, universal, atmos, layer, flux)
       if (fb == 0.0d0) exit
    end do
    root = b
!write(*,*) 'iter, b, fb', iter, b, fb
    if (iter == itmax) then
       write (*,*) 'ERROR!!! zbrent: Maximum number of interations exceeded'
       write (*,*) 'called from: ',msg
    end if

  end function zbrent

  !-----------------------------------------------------------------------
  subroutine quadratic (a, b, c, r1, r2)
    !
    ! !DESCRIPTION:
    ! Solve a quadratic equation for its two roots
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: a,b,c       ! Terms for quadratic equation
    real(8), intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    real(8) :: q                        ! Temporary term for quadratic solution
    !---------------------------------------------------------------------

    if (a == 0.0d0) then
       write (*,*) 'Quadratic solution error: a = ',a
    end if

    if (b >= 0.0d0) then
       q = -0.5d0 * (b + sqrt(b*b - 4.0d0*a*c))
    else
       q = -0.5d0 * (b - sqrt(b*b - 4.0d0*a*c))
    end if

    r1 = q / a
    if (q /= 0.0d0) then
       r2 = c / q
    else
       r2 = 1.e36
    end if

  end subroutine quadratic

  !-----------------------------------------------------------------------
  subroutine tridiag (a, b, c, r, u, n)
    !
    ! !DESCRIPTION:
    ! Solve a tridiagonal system of equations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: n        ! Number of soil layers
    real(8), intent(in)  :: a(n)     ! A vector for tridiagonal solution
    real(8), intent(in)  :: b(n)     ! B vector for tridiagonal solution
    real(8), intent(in)  :: c(n)     ! C vector for tridiagonal solution
    real(8), intent(in)  :: r(n)     ! R vector for tridiagonal solution
    real(8), intent(out) :: u(n)     ! U vector for tridiagonal solution
    !
    ! !LOCAL VARIABLES:
    real(8) :: gam(n)                ! Temporary calculation
    real(8) :: bet                   ! Temporary calculation
    integer :: j                     ! Soil layer index
    !---------------------------------------------------------------------

    ! Tridiagonal solution:
    !
    ! Solve for U given the set of equations F x U = R, where U is a vector
    ! of length N, R is a vector of length N, and F is an N x N tridiagonal
    ! matrix defined by the vectors A, B, C (each of length N). A(1) and
    ! C(N) are undefined and are not referenced by the subroutine.
    !
    !    | b(1) c(1)   0  ...                      |   | u(1)   |   | r(1)   |
    !    | a(2) b(2) c(2) ...                      |   | u(2)   |   | r(2)   |
    !    |                ...                      | x | ...    | = | ...    |
    !    |                ... a(n-1) b(n-1) c(n-1) |   | u(n-1) |   | r(n-1) |
    !    |                ...   0    a(n)   b(n)   |   | u(n)   |   | r(n)   |
    !

    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j)*gam(j)
       u(j) = (r(j) - a(j)*u(j-1)) / bet
    end do
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do

  end subroutine tridiag

end module mod_math_tools
