! Coded on September 5, 2015
module outerinnerproduct
implicit none
contains
!### outer
function outer(v1,v2) result(v3)
real(8),intent(in) :: v1(3), v2(3)
real(8) :: v3(3)
v3(1)=v1(2)*v2(3) - v1(3)*v2(2)
v3(2)=v1(3)*v2(1) - v1(1)*v2(3)
v3(3)=v1(1)*v2(2) - v1(2)*v2(1)
return
end function outer
!### inner
function inner(v1,v2) result(v3)
real(8),intent(in) :: v1(3), v2(3)
real(8) :: v3
v3=v1(1)*v2(1) + v1(2)*v2(2) +v1(3)*v2(3)
return
end function inner
!### inner_rc
function inner_rc(v1,v2) result(v3)
real(8),   intent(in) :: v1(3)
complex(8),intent(in) :: v2(3)
complex(8) :: v3
v3=v1(1)*v2(1) + v1(2)*v2(2) +v1(3)*v2(3)
return
end function inner_rc
end module outerinnerproduct
