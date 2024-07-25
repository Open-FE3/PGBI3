module horizontalresolution
implicit none

contains
!########################################## value
function value(x,y,g_meshpara)
use param_mesh
implicit none
type(meshpara),intent(in) :: g_meshpara
real(8),       intent(in) :: x,y
real(8) :: si,sb, x1,x2,y1,y2,sigma,A
real(8) :: value, robs, xyz(3,g_meshpara%nobs)
integer(4) :: i,nobs

!#[0]## set input
nobs = g_meshpara%nobs
xyz  = g_meshpara%xyz
x1 = g_meshpara%xbound(2)
x2 = g_meshpara%xbound(3)
y1 = g_meshpara%ybound(2)
y2 = g_meshpara%ybound(3)
si = g_meshpara%sizein
sb = g_meshpara%sizebo
A = g_meshpara%A ! [km]
sigma = g_meshpara%sigma

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. y .lt. y1 .or. y2 .lt. y) then
  value=sb
 else
!#[2]## value for focus area or near observatories
  value=si
  do i=1,nobs
   robs=dsqrt((xyz(1,i)-x)**2.d0 + (xyz(2,i)-y)**2.d0 )
   if ( robs .le. 3.*sigma ) then
    value = min(value, si - (si - A)*exp(-robs**2./2./sigma**2.) )
   else
    value = min(value,si)
   end if
  end do
 end if

return
end

!########################################## value
function value_3d(x,y,z,g_meshpara)
use param_mesh
implicit none
type(meshpara),intent(in) :: g_meshpara
real(8),       intent(in) :: x,y,z
real(8) :: si,sb, x1,x2,y1,y2,z1,z2,sigma,A
real(8) :: value_3d, robs, xyz(3,g_meshpara%nobs)
integer(4) :: i,nobs

!#[0]## set input
nobs = g_meshpara%nobs
xyz  = g_meshpara%xyz
x1 = g_meshpara%xbound(2)
x2 = g_meshpara%xbound(3)
y1 = g_meshpara%ybound(2)
y2 = g_meshpara%ybound(3)
z2 = g_meshpara%upzin
z1 = g_meshpara%downzin
si = g_meshpara%sizein3d
sb = g_meshpara%sizebo
A  = g_meshpara%A ! [km]
sigma = g_meshpara%sigma

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. &
   &  y .lt. y1 .or. y2 .lt. y .or. &
   &  z .lt. z1 .or. z2 .lt. z ) then
  value_3d=sb
 else
!#[2]## value for focus area or near observatories
  value_3d=si
  do i=1,nobs
   robs=dsqrt((xyz(1,i)-x)**2.d0 + (xyz(2,i)-y)**2.d0 + (xyz(3,i)-z)**2.d0)
   if ( robs .le. 3.*sigma ) then
    value_3d = min(value_3d, si - (si - A)*exp(-robs**2./2./sigma**2.) )
   else
    value_3d = min(value_3d,si)
   end if
  end do
 end if

return
end

!########################################### upgradegrd
subroutine updategrd(xgrd,nxmax,nx,x,dd)
implicit none
integer(4),intent(in) :: nxmax
real(8),intent(inout) :: xgrd(nxmax)
integer(4),intent(inout) :: nx
real(8),intent(in) :: dd,x
integer(4) :: k,kk
kk=0
!write(*,*) "input y=",x
 do k=1,nx
  if ( xgrd(k) .lt. x ) kk=k ! remember k
  if ( dabs(xgrd(k) - x) .lt. dd ) goto 10 ! do nothing
 end do
 xgrd(kk+2:nx+1)=xgrd(kk+1:nx) ! shift
 xgrd(kk+1)=x ! new x
 nx=nx+1
! write(*,*) "new y =",xgrd(kk+1)
 if (nxmax .lt. nx) then
  write(*,*) "GEGEGE nx is greater than nx! nxmax=",nxmax,"nx=",nx
  stop
 end if
10 continue
return
end

end module horizontalresolution
