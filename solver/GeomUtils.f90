! File: GeomUtils.f90
! Purpose: Encapsulate module 'GeomUtils', which provides some
!          trigonometric-conversions (between degrees and radians), for
! NOTE: This is a simplified version of the 'GeomUtil' module from Chapter4
!       (box-model).
module GeomUtils
  use NumericKinds
  implicit none
  public

  real(RK), parameter :: &
       PI = 3.141593_RK, &
       ONE_DEG_IN_RADS = PI / 180._RK

contains
  ! Convert degrees to radians
  real(RK) function deg2Rad( degrees )
    real(RK), intent(in) :: degrees
    deg2Rad = degrees * ONE_DEG_IN_RADS
  end function deg2Rad

  ! Sine of an angle given in degrees
  real(RK) function sinD( degrees )
    real(RK), intent(in) :: degrees
    sinD = sin( deg2Rad( degrees ) )
  end function sinD
end module GeomUtils
