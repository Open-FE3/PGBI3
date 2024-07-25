!# list are changed on Sep 14, 2016
!#coded on Jan. 17, 2016
!# to summarize the info about mesh
module mesh_type
implicit none

type mesh
 !# mesh info
 character(70) :: meshname
 integer(4) :: node ! # of node
 integer(4) :: ntet ! # of tetrahedron
 integer(4) :: ntri ! # of triangles
 integer(4) :: nlin ! # of lines
 integer(4) :: npoi ! # of points
 integer(4) :: icoordinateflag ! main coordinate flag, added on 2016.09.03
 !             icoordinateflag = 1 for xyz, 2 for lonlatalt, 3 for xyzspherical
 real(8)                               :: lonorigin,latorigin
 real(8),   allocatable,dimension(:,:) :: xyz
 real(8),   allocatable,dimension(:,:) :: lonlatalt    ! added on 2016.09.02
 real(8),   allocatable,dimension(:,:) :: xyzspherical ! added on 2016.09.02
 integer(4),allocatable,dimension(:,:) :: n4 ! 4 node-id for tetrahedrons
 integer(4),allocatable,dimension(:,:) :: n3 ! 3 node-id for triangles
 integer(4),allocatable,dimension(:,:) :: n2 ! 2 node-id for lines
 integer(4),allocatable,dimension(:,:) :: n1 ! 1 node-id for points
 integer(4),allocatable,dimension(:,:) :: n4flag ! n4flag(1:ntet,1:2) for tet
 integer(4),allocatable,dimension(:,:) :: n3flag ! n3flag(1:ntri,1:2) for tri
 integer(4),allocatable,dimension(:,:) :: n2flag ! n2flag(1:nlin,1:2) for tet
 integer(4),allocatable,dimension(:,:) :: n1flag ! n1flag(1:npoi,1:2) for tet
 !# model info
 integer(4)                          :: nmodel     ! 2018.11.14
 integer(4),allocatable,dimension(:) :: iele2model ! i-th ele corresponds to iele2model(i)-th model field
 real(8),   allocatable,dimension(:) :: cmodel     ! conductivity for i-th field
end type

type list_type
 integer(4) :: ne ! # of element
 integer(4),allocatable,dimension(:) :: list
end type

type grid_list_type
 integer(4) :: nx,ny,nz ! number of obex in x,y,z direction
 integer(4),allocatable,dimension(:,:) :: nelist
   ! nelist(i,1:3) 1:3 are grid index, i is the element index
 integer(4),allocatable,dimension(:)   :: ordered_ele
 integer(4),allocatable,dimension(:)   :: stack
 real(8),   allocatable,dimension(:)   :: xgrd,ygrd,zgrd
 ! xgrd(1:nx+1),ygrd(1:ny+1),zgrd(1:nz+1)
end type

contains
!############################################# allocate grid_list
subroutine allocate_grid_list(nx,ny,nz,ntet,glist)
type(grid_list_type),intent(inout) :: glist
integer(4),intent(in) :: nx,ny,nz,ntet

glist%nx=nx
glist%ny=ny
glist%nz=nz
allocate(glist%nelist(ntet,3))
allocate(glist%stack(0:nx*ny*nz))
allocate(glist%ordered_ele(ntet))
allocate(glist%xgrd(1:nx+1))
allocate(glist%ygrd(1:ny+1))
allocate(glist%zgrd(1:nz+1))
return
end

!############################################# gen grid for list
subroutine gengridforlist(xyzminmax,glist)
type(grid_list_type),intent(inout) :: glist
real(8),intent(in) :: xyzminmax(6)
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
integer(4) :: nx,ny,nz,i,j,k
xmin=xyzminmax(1)
xmax=xyzminmax(2)
ymin=xyzminmax(3)
ymax=xyzminmax(4)
zmin=xyzminmax(5)
zmax=xyzminmax(6)
nx=glist%nx
ny=glist%ny
nz=glist%nz

!#[1] generate x,y,z grid
do i=1,nz+1
 glist%zgrd(i)=(i-1)*(zmax-zmin)/float(nz) + zmin
! write(*,*) "i",i,"zgrd(i)=",glist%zgrd(i)
end do
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
!# lt to le on 2017.06.29
subroutine classifyelement2grd(h_mesh,glist)
type(mesh),             intent(in)     :: h_mesh
type(grid_list_type),   intent(inout)  :: glist
real(8),    allocatable,dimension(:)   :: order
integer(4), allocatable,dimension(:)   :: index
real(8),    allocatable,dimension(:,:) :: xyz
integer(4)  :: ii,i,j,k,ix,iy,iz,nx,ny,nz
real(8)     :: x3(3)
integer(4)  :: ipre,iflag

!#[0]## set
nx  = glist%nx
ny  = glist%ny
nz  = glist%nz
allocate(order(0:h_mesh%ntet),index(h_mesh%ntet))
allocate(xyz(3,h_mesh%node)) ! 2016.11.20
xyz = h_mesh%xyz ! 2016.11.20

!#[1]# fill glist%nelist and calculate order
do ii=1,h_mesh%ntet
 x3(:)=0.d0
 do j=1,4
  x3(1:3)=x3(1:3)+xyz(1:3,h_mesh%n4(ii,j))/4.d0
   ! x,y,z coordinate of center of gravity
 end do
 do k=1,glist%nz
  if (glist%zgrd(k) .lt. x3(3) .and. x3(3) .le. glist%zgrd(k+1)) then
    iz=k
  do j=1,glist%ny
  if (glist%ygrd(j) .lt. x3(2) .and. x3(2) .le. glist%ygrd(j+1)) then
    iy=j
   do i=1,glist%nx
   if (glist%xgrd(i) .lt. x3(1) .and. x3(1) .le. glist%xgrd(i+1)) then
    ix=i ; goto 100
   end if
   end do
  end if
  end do
 end if
 end do
 write(*,*) "GEGEGE no corresponding element list...in making nelist"
 write(*,*) "x3(1:3)=",x3(1:3)
 stop
 100 continue
 glist%nelist(ii,1:3)=(/ix,iy,iz/)
 order(ii)= ((nx*ny)*(iz-1) + nx*(iy-1) + ix )*1.d0
 index(ii)= ii
end do

!#[2]# calculate glist%ordered_ele
CALL SORT_INDEX(h_mesh%ntet,index,order(1:h_mesh%ntet))
glist%ordered_ele = index

!#[3]# calculate glist%stack
glist%stack(:)=0 ; ipre=0 !## ipre is # of accumulated elements
do ii = 0, nx*ny*nz-1
 10 continue
 if ( ipre .eq. h_mesh%ntet ) goto 20
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
glist%stack(ii:nx*ny*nz)=h_mesh%ntet

write(*,*) "### CLASSIFYELEMENT2GRD END!! ###"
!do ii=-400,0
!write(*,*) "ii=",ii,"glist%stack(nxyz+ii)=",glist%stack(nx*ny*nz+ii),"order(ntet+ii)=",order(h_mesh%ntet+ii)
!end do
!stop

return
end
!############################################# findelewithgrid
subroutine findelewithgrid(h_mesh,glist,x3,iele,a4)
type(mesh),intent(in) :: h_mesh
type(grid_list_type),intent(in) :: glist
real(8),intent(in) :: x3(3)
integer(4),intent(out) :: iele
real(8),intent(out) :: a4(4)
integer(4) :: nx,ny,nz,i,j,k,l,ix,iy,iz,error,iside,index,num
type(list_type) :: elist
elist%ne=0
allocate(elist%list(h_mesh%ntet))
nx=glist%nx
ny=glist%ny
nz=glist%nz
  !#[1]## search corresponding grid for x3
  do l=1,nz
   if (glist%zgrd(l) .lt. x3(3) .and. x3(3) .le. glist%zgrd(l+1)) then
    iz=l
  do k=1,ny
   if (glist%ygrd(k) .lt. x3(2) .and. x3(2) .le. glist%ygrd(k+1)) then
    iy=k
   do j=1,nx
    if (glist%xgrd(j) .lt. x3(1) .and. x3(1) .le. glist%xgrd(j+1)) then
    ix=j ; goto 101
    end if
   end do
   end if
  end do
  end if
  end do
  write(*,*) "GEGEGE no corresponding element list...for x3(1:3)"
  write(*,*) "x3(1:3)=",x3(1:3)
  stop
  101 continue

  !#[2]## construct sum of nelist for search x3(1:3)
  iside=0
  100 continue ! return to here if error .ne. 0
  elist%ne=0
  do l=max(iz-iside,1),min(iz+iside,nz)
   do k=max(iy-iside,1),min(iy+iside,ny)
    do j=max(ix-iside,1),min(ix+iside,nx)
    if ( iz - (iside-1) .le. l .and. l .le. iz + (iside-1) .and. &
    &    iy - (iside-1) .le. k .and. k .le. iy + (iside-1) .and. &
    &    ix - (iside-1) .le. j .and. j .le. ix + (iside-1)  ) goto 20
      index=(nx*ny)*(l-1) + nx*(k-1) + j
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
  end do
  if (  elist%ne .eq. 0 ) then
    iside=iside+1
     goto 100 ! iside should be increased
  end if
  !#[3]## search x3 in the nelist_sum
!  write(*,*) "ix,iy,iz=",ix,iy,iz
  call FINDELEMENT(x3,h_mesh,iele,a4,1,elist%ne,elist,error) ! see below
!
if (error .ne. 0 ) then ! if element was not , error = -1
! write(*,*) "iside=",iside, "to",iside+1
 iside=iside+1
 goto 100
end if

return
end

!---------------------------------------------  READMESH_TOTAL
! added on May 16, 2016
subroutine READMESH_TOTAL(h_mesh,infile)
implicit none
type(mesh),intent(inout) :: h_mesh
character(70),intent(in) :: infile

CALL MSHCOUNT1(infile,h_mesh)
CALL ALLOCATE_MESH(h_mesh)
CALL READMSH2(infile,h_mesh)

return
end
!---------------------------------------------  allocate_mesh
subroutine allocatemeshspherical(h_mesh)
implicit none
type(mesh),intent(inout) :: h_mesh

!# lonlatalt and xyzspherical
allocate (h_mesh%lonlatalt(3,h_mesh%node))
allocate (h_mesh%xyzspherical(3,h_mesh%node))

return
end
!---------------------------------------------  allocate_mesh
subroutine allocate_mesh(h_mesh)
implicit none

type(mesh),intent(inout) :: h_mesh

!# xyz
allocate (h_mesh%xyz(3,h_mesh%node))
!# node composing elements
allocate (h_mesh%n4(h_mesh%ntet,4))
allocate (h_mesh%n3(h_mesh%ntri,3))
allocate (h_mesh%n2(h_mesh%nlin,2))
allocate (h_mesh%n1(h_mesh%npoi,1))

!# flags
allocate (h_mesh%n4flag(h_mesh%ntet,2))
allocate (h_mesh%n3flag(h_mesh%ntri,2))
allocate (h_mesh%n2flag(h_mesh%nlin,2))
allocate (h_mesh%n1flag(h_mesh%npoi,2))

!# model
!allocate (h_mesh%iele2model(h_mesh%ntet)) ! commented out on 2018.11.14
allocate (h_mesh%cmodel(h_mesh%ntet))     ! 2018.11.14

write(*,*) "allocate mesh end!! ",h_mesh%meshname
end

!########################################### mshcount1
subroutine mshcount1(mshfile,h_mesh)
implicit none
type(mesh),intent(out) :: h_mesh
integer(4)             :: node,ntet,ntri,nlin,npoi
integer(4)             :: i,j,k,itype,ii,jj,kk,nele,n44(4)
!integer(4) :: etype
character(70),intent(in) :: mshfile
!# set the filename as meshname
h_mesh%meshname=mshfile

!# open the file
open(1,file=mshfile)
do i=1,4
read(1,*)
end do
read(1,*) node
!write(*,*) "node=",node
do i=1,node+2 ! include "$EndNodes", "$Elements" lines
read(1,*)
end do
read(1,*) nele ! total number of elements
npoi=0;nlin=0;ntri=0;ntet=0
do i=1,nele
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
npoi=npoi+1
else if ( itype .eq. 1) then ! Line element
nlin=nlin+1
else if ( itype .eq. 2) then   ! Triangle element
ntri=ntri+1
else if ( itype .eq. 4) then ! Tetrahedral element
ntet=ntet+1
end if
end do
close(1)
h_mesh%node=node
h_mesh%ntet=ntet
h_mesh%ntri=ntri
h_mesh%nlin=nlin
h_mesh%npoi=npoi
!write(*,*) "# of Point elements is",npoi
!write(*,*) "# of Line elements is",nlin
!write(*,*) "# of Triangle elements is",ntri
!write(*,*) "# of Tetrahedron elements is",ntet
write(*,*) "### COUNT ", mshfile(1:len_trim(mshfile))," END!! ###"
return
end subroutine mshcount1
!##########################################  readmsh2
subroutine readmsh2(mshfile,h_mesh)
implicit none
type(mesh),   intent(inout) :: h_mesh
character(70),intent(in)    :: mshfile
integer(4)                  :: itet,itri,ilin,ipoi,nele,inode
integer(4)                  :: i,j,k,ii,jj,kk,itype,nkk,ntet2,ifile
!integer(4) :: etype
integer(4),dimension(4) :: n44
open(1,file=mshfile)
!# [1] ### skip header
do i=1,4
read(1,*)
end do
!# [2] ### read node
read(1,*) inode
write(*,*) "# of nodes (node) =",inode
if ( inode .gt. h_mesh%node) goto 999
do i=1,h_mesh%node
read(1,*) j,(h_mesh%xyz(k,i),k=1,3) ! [km]
end do
read(1,*)
!# [3] ### read elements
read(1,*) ! skip the starting line, "$Elements"
read(1,*) nele
write(*,*)"# of elements (nele)=",nele
ipoi=0;ilin=0;itri=0;itet=0;nkk=0
do i=1,nele
! jj: pysical group, kk: volume id
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
ipoi=ipoi+1
h_mesh%n1    (ipoi,1)=n44(1)
h_mesh%n1flag(ipoi,1:2)=(/jj,kk/)
else if ( itype .eq. 1) then ! Line element
ilin=ilin+1
h_mesh%n2    (ilin,1:2)=n44(1:2)
h_mesh%n2flag(ilin,1:2)=(/jj,kk/)
else if ( itype .eq. 2 ) then ! Triangle element
itri=itri+1
h_mesh%n3    (itri,1:3)=n44(1:3)
h_mesh%n3flag(itri,1:2)=(/jj,kk/)
if (kk .eq. 1) nkk=nkk+1 ! count triangle group 1
else if ( itype .eq. 4) then ! Tetrahedral element
itet=itet+1
h_mesh%n4    (itet,1:4)=n44(1:4)
h_mesh%n4flag(itet,1:2)=(/jj,kk/)
end if
end do
!------------------------------------------------
write(*,*) "# of Point elements is",      h_mesh%npoi
write(*,*) "# of Line elements is",       h_mesh%nlin
write(*,*) "# of Triangle elements is",   h_mesh%ntri
write(*,*) "# of traiangle group1 is",    nkk
write(*,*) "# of Tetrahedron elements is",h_mesh%ntet
read(1,*) ! "END ELEMENT"
read(1,*,end=99) ! file end or $ElementData
ifile=1
call readcmodel(ifile,h_mesh)
99 continue
close(1)
write(*,*) "### READ ", mshfile(1:len_trim(mshfile))," END!! ###"
return
999 write(*,*) "GEGEGE node .ne. inode, node=",h_mesh%node,"inode=",inode
stop
end subroutine readmsh2

!############################################## read conductivity model
subroutine readcmodel(ifile,h_mesh)
implicit none
integer(4),intent(in) :: ifile
type(mesh),intent(inout) :: h_mesh
integer(4) :: ntet2, ii,i
character(80) a
do i=1,7
 read(ifile,'(a)') a
! write(*,*) a
end do
read(ifile,*) ntet2
write(*,*) ntet2
if ( ntet2 .ne. h_mesh%ntet ) goto 998
do i=1,h_mesh%ntet
 read(ifile,*) ii, h_mesh%cmodel(i)
end do
return
!
998 write(*,*) "GEGEGE ntet2 .ne. ntet, ntet=",h_mesh%ntet,"ntet2=",ntet2
stop
end subroutine readcmodel
!##
!
!######################################## OUTNODECVECTOR
!# coded on Nov. 27, 2015
!# iflag=1 : real, 2: imag
subroutine OUTNODECVECTOR(filename,vec, h_mesh,iflag)
implicit none
character(70),intent(in) :: filename
type(mesh),   intent(in) :: h_mesh
complex(8),   intent(in) :: vec(h_mesh%node)
integer(4),   intent(in) :: iflag
integer(4) :: i
real(8)    :: rvec(h_mesh%node)

!#[1]## check iflag
if (iflag .ne. 1 .and. iflag .ne. 2 )then
 write(*,*) "GEGEGE! iflag should be 1 (real) or 2(imag), iflag=",iflag
 stop
end if

!#[2]## generate real vector
 if (iflag .ne. 1 ) then ! real
  do i=1,h_mesh%node
   rvec(i)=dreal(vec(i))
  end do
 else                    ! imaginary
  do i=1,h_mesh%node
   rvec(i)=dimag(vec(i))
  end do
 end if

!#[3]## out real vector
CALL OUTNODERVECTOR(filename,rvec,h_mesh)

return
end
!######################################## OUTNODECVECTOR
subroutine OUTNODERSCALAR(filename,scalar,h_mesh)
implicit none
character(70),intent(in) :: filename
type(mesh),   intent(in) :: h_mesh
real(8),   intent(in) :: scalar(h_mesh%node)
integer(4) :: ifile,j

!#[1]# open file
ifile=1
open(ifile,file=filename)

!#[2]# out mesh info
CALL MESHOUT(ifile,h_mesh)

!#[3]## create vector fields
write(ifile,'(a)') "$NodeData"
write(ifile,'(a)') "1"
write(ifile,'(a)') '"A vector view"'
write(ifile,'(a)') "1"
write(ifile,'(a)') "0.0"
write(ifile,'(a)') "3"
write(ifile,'(a)') "0"
write(ifile,'(a)') "1" ! means only one (scalar) value is assigned to element
write(ifile,*) h_mesh%node
 write(1,'(i10,g15.7)') (j,scalar(j),j=1,h_mesh%node)
write(1,'(a)') "$EndNodeData"

!#[4]## cloase file
close(ifile)

write(*,*)"### OUTPUT TO ",filename(1:len_trim(filename))," END! ###"
return
end

!######################################## OUTNODERVECTOR
subroutine OUTNODERVECTOR(filename,rvec, h_mesh)
implicit none
character(70),intent(in) :: filename
type(mesh),   intent(in) :: h_mesh
real(8),      intent(in) :: rvec(3,h_mesh%node)
integer(4) :: ifile,i,j

!#[1]# open file
ifile=1
open(ifile,file=filename)

!#[2]# out mesh info
CALL MESHOUT(ifile,h_mesh)

!#[3]## create vector fields
write(ifile,'(a)') "$NodeData"
write(ifile,'(a)') "1"
write(ifile,'(a)') '"A vector view"'
write(ifile,'(a)') "1"
write(ifile,'(a)') "0.0"
write(ifile,'(a)') "3"
write(ifile,'(a)') "0"
write(ifile,'(a)') "3" ! means three vector is assigned to element
write(ifile,*) h_mesh%node
 write(1,'(i10,3g15.7)') ((j,rvec(i,j),i=1,3),j=1,h_mesh%node)
write(1,'(a)') "$EndNodeData"

!#[4]## cloase file
close(ifile)

write(*,*)"### OUTPUT TO ",filename(1:len_trim(filename))," END! ###"
return
end

!######################################## MESHOUT
subroutine MESHOUT(ifile,h_mesh)
implicit none
integer(4),intent(in) :: ifile
type(mesh),intent(in) :: h_mesh
integer(4) :: i,ishift

!#[0] Header
CALL HEADEROUT(ifile)

!#[1] Nodes
CALL NODEOUT(ifile,h_mesh)

!#[2] Elements
CALL ELEMENTOUT(ifile,h_mesh)

return
end
!########################################## HEADEROUT
subroutine HEADEROUT(ifile)
integer(4),intent(in) :: ifile
write(ifile,'(a)')"$MeshFormat"
write(ifile,'(a)')"2.2 0 8"
write(ifile,'(a)')"$EndMeshFormat"
return
end
!########################################## NODEOUT
subroutine NODEOUT(ifile,h_mesh)
type(mesh),intent(in) :: h_mesh
integer(4),intent(in) :: ifile
integer(4) :: i
real(8) :: xyz(3,h_mesh%node) ! added on 2016.09.03
write(ifile,'(a)')"$Nodes"
write(ifile,*) h_mesh%node

!# set coordinate type, added on 2016.09.03
write(*,*) "coordinate type 1 for cartesian (surface origin)"
write(*,*) "                2 for lon,lat,altitude"
write(*,*) "                3 for spherical cartesian (earth center origin)"
write(*,*) "output coordinate type is", h_mesh%icoordinateflag
xyz(1:3,:)=h_mesh%xyz(1:3,:) ! default is normal xyz
if (h_mesh%icoordinateflag .eq. 1) xyz(1:3,:)=h_mesh%xyz(1:3,:)
if (h_mesh%icoordinateflag .eq. 2) xyz(1:3,:)=h_mesh%lonlatalt(1:3,:)
if (h_mesh%icoordinateflag .eq. 3) xyz(1:3,:)=h_mesh%xyzspherical(1:3,:)

do i=1,h_mesh%node
 write(ifile,*) i, xyz(1:3,i)
end do
write(ifile,'(a)')"$EndNodes"
return
end
!########################################## ELEMENTOUT
subroutine ELEMENTOUT(ifile,h_mesh)
type(mesh),intent(in) :: h_mesh
integer(4),intent(in) :: ifile
integer(4) :: i,ishift
write(ifile,'(a)')"$Elements"
write(ifile,*) h_mesh%npoi + h_mesh%nlin + h_mesh%ntri + h_mesh%ntet
ishift=0
do i=1,h_mesh%npoi
 write(ifile,*) i + ishift," 15 2 ",h_mesh%n1flag(i,1:2),h_mesh%n1(i,1)
end do
ishift=h_mesh%npoi
do i=1,h_mesh%nlin
 write(ifile,*) i + ishift," 1 2 ",h_mesh%n2flag(i,1:2),h_mesh%n2(i,1:2)
end do
ishift=ishift + h_mesh%nlin
do i=1,h_mesh%ntri
 write(ifile,*) i + ishift," 2 2 ",h_mesh%n3flag(i,1:2),h_mesh%n3(i,1:3)
end do
ishift=ishift + h_mesh%ntri
do i=1,h_mesh%ntet
 write(ifile,*) i + ishift," 4 2 ",h_mesh%n4flag(i,1:2),h_mesh%n4(i,1:4)
end do
write(ifile,'(a)')"$EndElements"

return
end

!######################################## FINDELEMENT
! istart and iend are added on May 17, 2016
! This subroutine finds corresponding elements that includes given coord
subroutine FINDELEMENT(x3,h_mesh,iele,a,istart,iend,nelist,error)
!use mesh_type
use fem_util ! for volume, see m_fem_util.f90
implicit none
type(mesh),intent(in)  :: h_mesh
integer(4),optional,intent(in) :: istart,iend ! added on May 17, 2016
type(list_type),optional,intent(in) :: nelist ! added on May 17, 2016
real(8),   intent(in)  :: x3(3) ! [km]
integer(4),intent(out) :: iele
real(8)   ,intent(out) :: a(4)
integer(4) :: i, j,n(4),is,ie,ii
integer(4),intent(out) :: error
real(8) :: elm_xyz(3,4),x3_center(3)
is=1;ie=h_mesh%ntet            ! added on May 17, 2016
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
  n(1:4)=h_mesh%n4(ii,1:4) ! note that n4(i,1) indicates element group
  x3_center(1:3)=0.d0
  do j=1,4
   elm_xyz(1:3,j)=h_mesh%xyz(1:3,n(j))
   x3_center(1:3)=x3_center(1:3) + elm_xyz(1:3,j)/4.d0
  end do
!  write(*,*) "i=",i,"x3_center(1:3)=",x3_center(1:3),"is=",is,"ie=",ie
  CALL nodebasisfun(elm_xyz,x3,a) ! see m_fem_util.f90

 !#[4]## check weather the point is included in i-th element
 if ( a(1) .ge. 0.d0 .and. a(2) .ge. 0.d0   .and. &
    & a(3) .ge. 0.d0 .and. a(4) .ge. 0.d0 ) then
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
!write(*,*) "GEGEGE! No corresponding element was found! return"
!write(*,*) "x,y,z=",x3(1:3)
!write(*,*) "is=", is
!write(*,*) "ie=", ie
!stop
end

!###########################################  function etype
function etype(itype)
implicit none
integer(4),intent(in) :: itype
integer(4) :: etype
etype=0
if (itype .eq. 15) etype=1! one point node
if (itype .eq. 1 ) etype=2 ! line
if (itype .eq. 2 ) etype=3 ! triangle
if (itype .eq. 4 ) etype=4 ! tetrahedron
return
end function etype
!

end module mesh_type


