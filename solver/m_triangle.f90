! coded on 2016.11.20
module triangle
use mesh_type
implicit none



contains

!############################################# allocate grid_list
subroutine allocate_2Dgrid_list(nx,ny,ntri,glist)
type(grid_list_type),intent(inout) :: glist
integer(4),intent(in) :: nx,ny,ntri

glist%nx=nx
glist%ny=ny
allocate(glist%nelist(ntri,2))
allocate(glist%stack(0:nx*ny))
allocate(glist%ordered_ele(ntri))
allocate(glist%xgrd(1:nx+1))
allocate(glist%ygrd(1:ny+1))
return
end

!############################################# gen grid for list
subroutine gen2Dgridforlist(xyzminmax,glist)
type(grid_list_type),intent(inout) :: glist
real(8),intent(in) :: xyzminmax(6)
real(8) :: xmin,xmax,ymin,ymax
integer(4) :: nx,ny,i,j
xmin=xyzminmax(1)
xmax=xyzminmax(2)
ymin=xyzminmax(3)
ymax=xyzminmax(4)
nx=glist%nx
ny=glist%ny

!#[1] generate x,y,z grid
do i=1,ny+1
 glist%ygrd(i)=(i-1)*(ymax-ymin)/float(ny) + ymin
! write(*,*) "i",i,"ygrd(i)=",glist%ygrd(i)
end do
do i=1,nx+1
 glist%xgrd(i)=(i-1)*(xmax-xmin)/float(nx) + xmin
! write(*,*) "i",i,"xgrd(i)=",glist%xgrd(i)
end do

return
end

!############################################# classifyelement2grd
subroutine classifytri2grd(h_mesh,glist)
type(mesh),            intent(in)     :: h_mesh
type(grid_list_type),  intent(inout)  :: glist
real(8)                               :: x2(2)
integer(4)                            :: ii,i,j,ix,iy,nx,ny
real(8),   allocatable,dimension(:)   :: order
integer(4),allocatable,dimension(:)   :: index
real(8),   allocatable,dimension(:,:) :: xyz
integer(4)                            :: ipre,iflag

!#[0]## set
nx=glist%nx
ny=glist%ny
allocate(order(0:h_mesh%ntri),index(h_mesh%ntri))
allocate(xyz(2,h_mesh%node)) ! 2016.11.20
xyz(:,:) = h_mesh%xyz(1:2,:) ! 2016.11.20

!#[1]# fill glist%nelist and calculate order
do ii=1,h_mesh%ntri
 x2(:)=0.d0
 do j=1,3
  x2(1:2)=x2(1:2)+xyz(1:2,h_mesh%n3(ii,j))/3.d0
   ! x,y,z coordinate of center of gravity
 end do
  do j=1,glist%ny
  if (glist%ygrd(j) .lt. x2(2) .and. x2(2) .lt. glist%ygrd(j+1)) then
    iy=j
   do i=1,glist%nx
   if (glist%xgrd(i) .lt. x2(1) .and. x2(1) .lt. glist%xgrd(i+1)) then
    ix=i ; goto 100
   end if
   end do
  end if
  end do
 write(*,*) "GEGEGE no corresponding element list...in making nelist"
 write(*,*) "x2(1:2)=",x2(1:2)
 stop
 100 continue
 glist%nelist(ii,1:2)=(/ix,iy/)
 order(ii)= ( nx*(iy-1) + ix )*1.d0
 index(ii)= ii
end do

!#[2]# calculate glist%ordered_ele
CALL SORT_INDEX(h_mesh%ntri,index,order(1:h_mesh%ntri))
glist%ordered_ele = index

!#[3]# calculate glist%stack
glist%stack(:)=0 ; ipre=0 !## ipre is # of accumulated elements
do ii = 0, nx*ny-1
 10 continue
 if ( ipre .eq. h_mesh%ntri ) goto 20
 if ( order(ipre) .lt. order(ipre+1) ) then
  glist%stack(ii)=ipre
  if ( ii + 1 .eq. int(order(ipre+1)) ) ipre = ipre + 1 ! if step is 1, add 1
 else                     !## when ( order(ipre) .eq. order(ipre+1) )
  ipre=ipre+1
  goto 10
 end if
! write(*,*) "ii=",ii,"glist%stack(ii)",glist%stack(ii)
end do
20 continue
glist%stack(ii:nx*ny)=h_mesh%ntri

!do ii=-400,0
!write(*,*) "ii=",ii,"glist%stack(nxyz+ii)=",glist%stack(nx*ny*nz+ii),"order(ntet+ii)=",order(h_mesh%ntet+ii)
!end do
!stop

return
end

!############################################# findtriwithgrid
!# based on findelewith grid
subroutine findtriwithgrid(h_mesh,glist,x2,iele,a3)
type(mesh),intent(in) :: h_mesh
type(grid_list_type),intent(in) :: glist
real(8),    intent(in) :: x2(2)
integer(4),intent(out) :: iele
real(8),   intent(out) :: a3(3)
integer(4) :: nx,ny,i,j,k,ix,iy,error,iside,index,num
type(list_type) :: elist
elist%ne=0
allocate(elist%list(h_mesh%ntri))
nx=glist%nx
ny=glist%ny
  !#[1]## search corresponding grid for x3
  do k=1,ny
   if (glist%ygrd(k) .lt. x2(2) .and. x2(2) .le. glist%ygrd(k+1)) then
    iy=k
   do j=1,nx
    if (glist%xgrd(j) .lt. x2(1) .and. x2(1) .le. glist%xgrd(j+1)) then
    ix=j ; goto 101
    end if
   end do
   end if
  end do
  write(*,*) "GEGEGE no corresponding element list...for x2(1:2)"
  write(*,*) "x2(1:2)=",x2(1:2)
  stop
  101 continue

  !#[2]## construct sum of nelist for search x3(1:3)
  iside=0
  100 continue ! return to here if error .ne. 0
  elist%ne=0
   do k=max(iy-iside,1),min(iy+iside,ny)
    do j=max(ix-iside,1),min(ix+iside,nx)
    if ( iy - (iside-1) .le. k .and. k .le. iy + (iside-1) .and. &
    &    ix - (iside-1) .le. j .and. j .le. ix + (iside-1)  ) goto 20
      index=nx*(k-1) + j
	num=glist%stack(index) - glist%stack(index-1)
!	write(*,*) "l,k,j=",l,k,j,"iside=",iside
!	write(*,*) "glist%stack(index-1, index)=",glist%stack(index-1:index)
!	write(*,*) "num=",num
	do i=1,num
       elist%list(elist%ne + i)=glist%ordered_ele(glist%stack(index-1)+i)
!	 write(*,*) "i=",i,"ordered_ele(glist%stack(index-1)+i)=",glist%ordered_ele(glist%stack(index-1)+i)
	end do
      elist%ne=elist%ne + num
    20 continue
    end do
   end do
  if (  elist%ne .eq. 0 ) then
    iside=iside+1
     goto 100 ! iside should be increased
  end if
  !#[3]## search x3 in the nelist_sum
!  write(*,*) "ix,iy,iz=",ix,iy,iz
  call FINDTRI(x2,h_mesh,iele,a3,1,elist%ne,elist,error) ! see below
!
if (error .ne. 0 ) then ! if element was not , error = -1
! write(*,*) "iside=",iside, "to",iside+1
 iside=iside+1
 goto 100
end if

return
end

!#####################################
! Coded on 2016.11.20 based on FINDELEMENT in m_mesh_type.f90
subroutine findtri(x2,h_mesh,iele,a,istart,iend,nelist,error)
!use mesh_type
use fem_util ! for volume, see m_fem_util.f90
implicit none
type(mesh),intent(in)  :: h_mesh
integer(4),optional,intent(in) :: istart,iend ! added on May 17, 2016
type(list_type),optional,intent(in) :: nelist ! added on May 17, 2016
real(8),   intent(in)  :: x2(2) ! [km]
integer(4),intent(out) :: iele
real(8)   ,intent(out) :: a(3)
integer(4) :: i, j,n(3),is,ie,ii
integer(4),intent(out) :: error
real(8) :: elm_xyz(2,3),a0
real(8),dimension(2) :: x12,x13,x23,v1,v2,v3
is=1;ie=h_mesh%ntri            ! added on May 17, 2016
if (present(nelist)) ie=nelist%ne
if (present(istart)) is=istart ! added on May 17, 2016
if (present(iend)  ) ie=iend   ! added on May 17, 2016
!
error=0

do i=is,iend ! modified on May 17, 2016
!if ( n4g(i,1) .eq. 1 ) then ! only elements under surface
ii=i
if (present(nelist) ) ii=nelist%list(i)

 !#[1]## generate n(1:4)
  n(1:3)=h_mesh%n3(ii,1:3) ! note that n4(i,1) indicates element group
  do j=1,3
   elm_xyz(1:2,j)=h_mesh%xyz(1:2,n(j))
  end do
 !# calculate a(1:3)
   x12(1:2) = elm_xyz(1:2,2) - elm_xyz(1:2,1) ! [km]
   x13(1:2) = elm_xyz(1:2,3) - elm_xyz(1:2,1) ! [km]
   x23(1:2) = elm_xyz(1:2,3) - elm_xyz(1:2,2) ! [km]
   v1 = x2(1:2) - elm_xyz(1:2,1)
   v2 = x2(1:2) - elm_xyz(1:2,2)
   v3 = x2(1:2) - elm_xyz(1:2,3)
   a0 =( x13(1)*x12(2) - x13(2)*x12(1))
   a(2)=( x13(1)* v1(2) - x13(2)* v1(1))/a0
   a(3)=(  v1(1)*x12(2) -  v1(2)*x12(1))/a0
   a(1)=(  v2(1)*x23(2) -  v2(2)*x23(1))/a0

 !#[4]## check weather the point is included in i-th element
 if ( a(1) .ge. 0.d0 .and. a(2) .ge. 0.d0 .and. a(3) .ge. 0.d0 ) then
  goto 100
 end if

end do
goto 999 !

!# the corresponding element was chosen
100 continue
iele=ii
return
!
999 continue
error=-1 ! not found
return
end

end module triangle
