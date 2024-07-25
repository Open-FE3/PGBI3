! coded on 2016.09.03
module constants
implicit none

! earthrad, pi, d2r, r2d
real(8),parameter :: earthrad=6371.d0, pi=4.d0*datan(1.d0),d2r=pi/180., r2d=180./pi

! mu
real(8),parameter :: dmu=4.d0*pi*1.d-7

end module constants
