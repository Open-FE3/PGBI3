module matrix
implicit none

type real_crs_matrix ! compressed row matrix
 integer(4) :: nrow
 integer(4) :: ntot
 integer(4) :: ncolm  ! 2017.11.02
 integer(4),dimension(:),allocatable :: item !item(1:ntot) :colmun id for every entry for
 integer(4),dimension(:),allocatable :: stack!stack(0:nrow) accumulated # of row members
 real(8),dimension(:),allocatable :: val  !val(1:ntot) : values
end type

contains
!#############################################  allocate_real_crs_with_steady_ncolm
!# coded on May 20, 2016
!# when the # of non zero members in every colmun
!# set mcrs%nrow, mcrs%ntot, and allocate mcrs%stack, mcrs%item, mcrs%val
!# and set mcrs%stack
subroutine allocate_real_crs_with_steady_ncolm(mcrs,nrow,ncolm)
implicit none
integer(4),intent(in) :: nrow,ncolm
type(real_crs_matrix),intent(inout) :: mcrs
integer(4) :: i
 mcrs%nrow  = nrow ! number of observation points
 mcrs%ntot  = nrow * ncolm ! for lines of tetrahedral mesh
 mcrs%ncolm = ncolm  ! 2017.11.02
 allocate(mcrs%stack(0:mcrs%nrow))
 allocate(mcrs%item(   mcrs%ntot))
 allocate(mcrs%val(    mcrs%ntot))
 mcrs%stack(0)=0
 do i=1,nrow
  mcrs%stack(i)=i*ncolm
 end do
return
end
!#############################################  mul_matcrs_rv
subroutine mul_matcrs_rv(crsmat,rv,doftot,vout)
implicit none
type(real_crs_matrix),intent(in)  :: crsmat
integer(4),           intent(in)  :: doftot
real(8),              intent(in)  :: rv(doftot)
real(8),              intent(out) :: vout(crsmat%nrow)
integer(4) :: irow,jtot
integer(4) :: chunk=1000 ! 2017.06.27
vout(:)=0.d0

!$OMP PARALLEL SHARED(crsmat,rv,vout,CHUNK) PRIVATE(irow,jtot) ! 2017.06.27
!$OMP DO SCHEDULE (STATIC,chunk)
do irow=1,crsmat%nrow
! write(*,*) "irow=",irow,"crsmat%nrow"
 do jtot=crsmat%stack(irow-1)+1,crsmat%stack(irow)
!  write(*,*) "jtot=",jtot,"crsmat%item(jtot)=",crsmat%item(jtot),"doftot=",doftot
  vout(irow)=vout(irow) + crsmat%val(jtot)*rv(crsmat%item(jtot))
 end do
end do
!$OMP END PARALLEL

return
end

!#############################################  mul_matcrs_rv
subroutine mul_matcrs_cv(crsmat,cv,doftot,vout)
implicit none
type(real_crs_matrix),intent(in) :: crsmat
integer(4),intent(in)   :: doftot
complex(8),  intent(in) :: cv(doftot)
complex(8),  intent(out):: vout(crsmat%nrow)
integer(4) :: irow,jtot
!#[1] decompose cv to rcv and icv
vout(:)=(0.d0,0.d0)
do irow=1,crsmat%nrow
! write(*,*) "irow=",irow,"crsmat%nrow=",crsmat%nrow
 do jtot=crsmat%stack(irow-1)+1,crsmat%stack(irow)
! write(*,*) "jtot=",jtot,"crsmat%item(jtot)=",crsmat%item(jtot),"cv(crsmat%item(jtot))=",cv(crsmat%item(jtot))
  vout(irow)=vout(irow) + crsmat%val(jtot)*cv(crsmat%item(jtot))
 end do
end do
return
end
!#############################################  mul_atb
subroutine  mul_atb( a, b, c,  aj_leng, bi_leng, cj_leng  )
!    multiply matrix 'c' = 'a't * 'b'
!      c(aj_leng,cj_leng) = a(bi_leng,aj_leng) * b(bi_leng,cj_leng)
implicit    none
integer(4) , intent(in)  ::  aj_leng !  2-nd element length of array 'a'
integer(4) , intent(in)  ::  bi_leng  ! 1-st element length of array 'b'
integer(4) , intent(in)  ::  cj_leng  !  2-nd element length of array 'c'
real (8), dimension(bi_leng, aj_leng), intent(in)  ::  a! input matrix 'a'
real (8), dimension(bi_leng, cj_leng), intent(in)  ::  b !   input matrix 'b'
real (8), dimension(aj_leng, cj_leng), intent(out) ::  c   !  input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,aj_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
  do k=1,bi_leng
  c(i,j) = c(i,j) + a(k,i)*b(k,j)
enddo
enddo
enddo
return
end subroutine  mul_atb!


!###############################################  mul_ab
subroutine  mul_ab ( a, b, c,  ai_leng, bi_leng, cj_leng  )
!    multiply matrix 'c' = 'a' * 'b'
!      c(ai_leng,cj_leng) = a(ai_leng,bi_leng) * b(bi_leng,cj_leng)
implicit    none
integer(4), intent(in)  ::  ai_leng ! 1-st element length of array 'a'
integer(4), intent(in)  ::  bi_leng ! 1-st element length of array 'b'
integer(4), intent(in)  ::  cj_leng !  2-nd element length of array 'c'
real (8), dimension(ai_leng, bi_leng), intent(in)  ::  a !   input matrix 'a'
real (8), dimension(bi_leng, cj_leng), intent(in)  ::  b  !  input matrix 'b'
real (8), dimension(ai_leng, cj_leng), intent(out) ::  c ! input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,ai_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
    do k=1,bi_leng
      c(i,j) = c(i,j) + a(i,k)*b(k,j)
    enddo
  enddo
enddo
return
end subroutine  mul_ab
!############################################### mul_abt
subroutine  mul_abt( a, b, c, ai_leng, bj_leng, cj_leng  )
!    multiply matrix 'c' = 'a' * 'b't
!      c(ai_leng,cj_leng) = a(ai_leng,bj_leng) * b(cj_leng,bj_leng)
implicit    none
integer(4), intent(in)  ::  ai_leng !  1-st element length of array 'a'
integer(4), intent(in)  ::  bj_leng !  2-nd element length of array 'b'
integer(4), intent(in)  ::  cj_leng !  2-nd element length of array 'c'
real   (8), dimension(ai_leng, bj_leng), intent(in)  ::  a!   input matrix 'a'
real   (8), dimension(cj_leng, bj_leng), intent(in)  ::  b !  input matrix 'b'
real   (8), dimension(ai_leng, cj_leng), intent(out) ::  c!  input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,ai_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
    do k=1,bj_leng
      c(i,j) = c(i,j) + a(i,k)*b(j,k)
    enddo
  enddo
enddo
return
end subroutine  mul_abt

end module matrix
