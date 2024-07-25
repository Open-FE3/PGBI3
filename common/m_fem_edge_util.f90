! coded on March 3, 2016
module fem_edge_util
implicit none

integer(4),dimension(3,2) :: lm ! added on 2016.09.08
data lm(1,1:2) /1,2/
data lm(2,1:2) /2,3/
data lm(3,1:2) /3,1/

integer(4),dimension(6,2) :: kl ! edge definition
data kl(1,1:2) /1,2/
data kl(2,1:2) /2,3/
data kl(3,1:2) /1,3/
data kl(4,1:2) /1,4/
data kl(5,1:2) /2,4/
data kl(6,1:2) /3,4/

integer(4),dimension(4,3) :: lmn ! face definition
data lmn(1,1:3) /2,3,4/
data lmn(2,1:3) /1,4,3/
data lmn(3,1:3) /1,2,4/
data lmn(4,1:3) /1,3,2/

integer(4),dimension(4,6) :: C ! rotation matrix
data C(1,1:6) / 0,  1,  0,  0, -1,  1/
data C(2,1:6) / 0,  0, -1,  1,  0, -1/
data C(3,1:6) / 1,  0,  0, -1,  1,  0/
data C(4,1:6) /-1, -1,  1,  0,  0,  0/
contains

!######################################## FACEBASISFUN
subroutine FACEBASISFUN(elm_xyz,x3,wfe,v)
use fem_util
use outerinnerproduct
use matrix
implicit none
real(8),intent(in)  :: elm_xyz(3,4)
real(8),intent(in)  :: x3(3)
real(8),intent(out) :: wfe(3,6), v
real(8) :: wf(3,4)
real(8) :: gn(3,4), lambda(4)
integer(4) :: i,l,m,n

!#[1]## vector interpolation function (scheme 2)
 call gradnodebasisfun(elm_xyz,gn,v) !gn=grad N and v is volume,    (m_fem_util.f90)
 call nodebasisfun(elm_xyz,x3,lambda)!lambda is node basis function (m_fem_util.f90)

!#[2]## calculate face basis function
!# w_f=2(lambda_l*(gn_m \times gn_n)&
!#     + lambda_n*(gn_l \times gn_m)&
!#     + lambda_m*(gn_n \times gn_l))
 do i=1,4
  l=lmn(i,1) ; m=lmn(i,2) ; n=lmn(i,3)
  wf(1:3,i)=2.d0*(&
&   lambda(l)*outer(gn(1:3,m),gn(1:3,n)) &
& + lambda(n)*outer(gn(1:3,l),gn(1:3,m)) &
& + lambda(m)*outer(gn(1:3,n),gn(1:3,l)) )
 end do

!#[3]## calculate wf*[C]
 call mul_ab(wf,C*1.d0,wfe,3,4,6)
 
return
end subroutine FACEBASISFUN

!######################################## EDGEBASISFUN
subroutine EDGEBASISFUN(elm_xyz,x3,w,len,v)
use fem_util
implicit none
real(8),intent(in)  :: elm_xyz(3,4)
real(8),intent(in)  :: x3(3)
real(8),intent(out) :: w(3,6), len(6), v
real(8) :: gn(3,4), lambda(4)
integer(4) :: i,k,l

!#[2]## vector interpolation function (scheme 2)
 call gradnodebasisfun(elm_xyz,gn,v) !gn=grad N and v is volume,    (m_fem_util.f90)
 call nodebasisfun(elm_xyz,x3,lambda)!lambda is node basis function (m_fem_util.f90)
 call edgelen(elm_xyz,len)
!#[3]## calculate vector basis function
 do i=1,6
  k=kl(i,1) ; l=kl(i,2)
  w(1:3,i)=lambda(k)*gn(1:3,l) - lambda(l)*gn(1:3,k)
 end do

return
end subroutine EDGEBASISFUN
!######################################## EDGELEN
subroutine EDGELEN(elm_xyz,len)
implicit none
real(8),intent(in)  :: elm_xyz(3,4) ! [km]
real(8),intent(out) :: len(6)       ! [km]
real(8) :: x3(3)
integer(4) :: i

do i=1,6
 x3(1:3)=elm_xyz(1:3,kl(i,2)) - elm_xyz(1:3,kl(i,1))
 len(i)=dsqrt(x3(1)**2.d0 + x3(2)**2.d0 + x3(3)**2.d0)
end do

return
end subroutine EDGELEN

!
end module fem_edge_util
