! Generated on March 3, 2016
module fem_util
implicit none
!### subroutine list ####
! volume : calculate volume of tetrahedron
! intv   : calculate volume integral of N*N in a tetrahedron
! nodebasisfun : node basis function for a point in terms of a tetrahedron.
! gradnodebasisfun  : calculate gradient N
contains

!#############################################
subroutine gradnodebasisfun(elm_xyz,gn,v)
use outerinnerproduct
implicit none
real(8),intent(in)  :: elm_xyz(3,4)! [km]
real(8),intent(out) :: gn(3,4)     ! 1/[km]
real(8),intent(out) :: v           ! [km]^3
real(8) :: xx(3,6)     ! [km]
integer(4) :: i

  call calxmn(elm_xyz,xx)
  call volume(xx(:,1), xx(:,5), -xx(:,4), v) ! volume of this element [km^3]
  if ( v .le. 0 ) goto 999
  ! gradient of nodal shape function, where outer vector, [1/km]
  gn(:,1)=1.d0/6.d0/v*outer(elm_xyz(:,4)-elm_xyz(:,2), elm_xyz(:,3)-elm_xyz(:,2))
  gn(:,2)=1.d0/6.d0/v*outer(elm_xyz(:,3)-elm_xyz(:,1), elm_xyz(:,4)-elm_xyz(:,1))
  gn(:,3)=1.d0/6.d0/v*outer(elm_xyz(:,4)-elm_xyz(:,1), elm_xyz(:,2)-elm_xyz(:,1))
  gn(:,4)=1.d0/6.d0/v*outer(elm_xyz(:,2)-elm_xyz(:,1), elm_xyz(:,3)-elm_xyz(:,1))

return
!########################################### when volume is less than zero
999 continue
write(*,*) "GEGEGE! volume is less than 0"
do i=1,4
 write(*,*)  "i=",i,"elm_xyz(i)=",elm_xyz(1:3,i)
end do
write(*,'(a,3g15.7)')  "xx(:,1)=",xx(:,1)
write(*,'(a,3g15.7)')  "xx(:,5)=", xx(:,5)
write(*,'(a,3g15.7)')  "-xx(:,4)=",-xx(:,4)
write(*,'(a,3g15.7)')  "v=",v
stop
end
!#############################################
subroutine calxmn(elm_xyz,xx)
implicit none
real(8),intent(in)  :: elm_xyz(3,4)
real(8),intent(out) :: xx(3,6)
  xx(:,1)=elm_xyz(:, 4) - elm_xyz(:, 3) ! [km]
  xx(:,2)=elm_xyz(:, 4) - elm_xyz(:, 1)
  xx(:,3)=elm_xyz(:, 2) - elm_xyz(:, 4)
  xx(:,4)=elm_xyz(:, 3) - elm_xyz(:, 2)
  xx(:,5)=elm_xyz(:, 1) - elm_xyz(:, 3)
  xx(:,6)=elm_xyz(:, 2) - elm_xyz(:, 1)
return
end
!#############################################
! if x3 is within the element, all of a are larger than 0.d0
subroutine nodebasisfun(elm_xyz,x3,a)
implicit none
real(8),intent(in) :: elm_xyz(3,4)
real(8),intent(in) :: x3(3)
real(8),intent(out) :: a(4)
integer(4) :: j
real(8),dimension(3) :: x12,x13,x14,x24,x23
real(8) :: a1,a2,a3,a4,a5,v(3,4)

!#[1]## vector from x3
 do j=1,4
  v(1:3,j)=elm_xyz(1:3,j) - x3(1:3) ! 4 vectors from point in concern
 end do

 x12(:)=elm_xyz(:,2) - elm_xyz(:,1)
 x13(:)=elm_xyz(:,3) - elm_xyz(:,1)
 x14(:)=elm_xyz(:,4) - elm_xyz(:,1)
 x24(:)=elm_xyz(:,4) - elm_xyz(:,2)
 x23(:)=elm_xyz(:,3) - elm_xyz(:,2)
 call volume(x12,x13,    x14, a5)
 call volume(x12,x13,-v(:,1), a4) ; a4=a4/a5
 call volume(x14,x12,-v(:,1), a3) ; a3=a3/a5
 call volume(x13,x14,-v(:,1), a2) ; a2=a2/a5
 call volume(x24,x23,-v(:,2), a1) ; a1=a1/a5
 a(1:4)=(/a1,a2,a3,a4/)

return
end
!#############################################
! integral of shape function
function intv(k, m, v)
implicit none
integer(4),intent(in) :: k, m
real(8),intent(in) :: v ![km^3]
real(8) :: intv
! here assume k = l
! intv =6v*(k!l!m!n!)/(k+l+m+n+3)!
if ( k .ne. m ) intv=v/20.d0 ![km^3]
if ( k .eq. m ) intv=v/10.d0 ![km^3]
return
end function intv

!#############################################
! assume closs produxt of x1 and x2 point on the side of direction of x3
! volume can be minus when (x1 times x2) cdot x3 < 0
subroutine volume(x1, x2, x3, v)
use outerinnerproduct
implicit none
real(8),dimension(3),intent(in) ::x1,x2,x3
real(8),intent(out) :: v
real(8),dimension(3) :: x12
v=1.d0/6.d0*inner(outer(x1,x2), x3)
return
end subroutine volume

!############################################# 4area
! Coded for 2017.05.18
subroutine area4(elm_xyz,a4)
use outerinnerproduct
implicit none
real(8),intent(in)  :: elm_xyz(3,4)
real(8),intent(out) :: a4(4)
real(8),dimension(3) :: x12,x13,x14,x23,x24,x34
x12 = elm_xyz(:,2) - elm_xyz(:,1)
x13 = elm_xyz(:,3) - elm_xyz(:,1)
x14 = elm_xyz(:,4) - elm_xyz(:,1)
x23 = elm_xyz(:,3) - elm_xyz(:,2)
x24 = elm_xyz(:,4) - elm_xyz(:,2)
a4(1) = dsqrt(inner(1./2.*outer(x23,x24),1./2.*outer(x23,x24)))
a4(2) = dsqrt(inner(1./2.*outer(x13,x14),1./2.*outer(x13,x14)))
a4(3) = dsqrt(inner(1./2.*outer(x12,x14),1./2.*outer(x12,x14)))
a4(4) = dsqrt(inner(1./2.*outer(x12,x13),1./2.*outer(x12,x13)))

return
end




end module fem_util
