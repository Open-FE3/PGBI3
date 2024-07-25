!# coded on 2017.10.25
module divroterror
use outerinnerproduct
use fem_edge_util
use param
use oceanFvxyz
use mesh_type
use fem_util
implicit none

type error
 ! rot error
 real(8)    :: x3_rot(3)
 real(8)    :: ediv_total
 real(8)    :: ediv_max
 integer(4) :: idev_div
 real(8)    :: dx3v3_rot(3,3)
 real(8)    :: gn_rot(3,4)
 real(8)    :: v3_rot(3,4)
 ! div error
 real(8)    :: x3_div(3)
 real(8)    :: erot_total
 real(8)    :: erot_max
 integer(4) :: idev_rot
 real(8)    :: dx3v3_div(3,3)
 real(8)    :: gn_div(3,4)
 real(8)    :: v3_div(3,4)
end type

contains

!########################################## prepoceanele2column
!# coded on 2017.11.02
subroutine outdivrotudistribution(g_param,em_mesh,h_ocean,it) ! 2017.11.08
implicit none
type(param_forward),   intent(in)     :: g_param
type(mesh),            intent(in)     :: em_mesh
type(ocean_data),      intent(in)     :: h_ocean
integer(4),            intent(in)     :: it
type(real_crs_matrix)                 :: oceanele2column
real(8),   allocatable,dimension(:)   :: diverr, roterr
real(8)                               :: ediv, erot
integer(4)                            :: iele
integer(4)                            :: i,j,k,iflag,jflag
integer(4)                            :: node, ntet, ntris
real(8),   allocatable,dimension(:,:) :: xyz
real(8),   allocatable,dimension(:,:) :: vxyz
integer(4),allocatable,dimension(:,:) :: n4
real(8),   allocatable,dimension(:,:) :: xy_tri
real(8)                               :: elm_xyz(3,4), v3(3,4)
real(8)                               :: denominator,v3abs
real(8)                               :: v3threshold = 1.0 ! [mm/sec]
character(70)                         :: outputfolder
character(70)                         :: outputfile
character(6)                          :: num

!#[0]## set
 outputfolder    = g_param%outvxyzfolder    ! 2018.11.14
 oceanele2column = h_ocean%oceanele2column
 node            = em_mesh%node
 ntet            = em_mesh%ntet
 allocate( xyz(3,node), n4(ntet,4), vxyz(3,node) )
 xyz             = em_mesh%xyz
 n4              = em_mesh%n4
 vxyz            = h_ocean%vxyz
 ntris           = h_ocean%ntris
 allocate( xy_tri(2,ntris) )
 xy_tri          = h_ocean%xy_tri

!#[1]## obtain diverr for every surface triangles
 allocate( diverr(ntris), roterr(ntris) )
 diverr = -9999.0
 roterr = -9999.0 ! 2017.11.08
 do i=1,ntris
  do j=  oceanele2column%stack(i-1)+1,oceanele2column%stack(i)
   iele = oceanele2column%item(j)
   do k=1,4
    v3(     1:3,k) = vxyz(1:3,n4(iele,k))
    elm_xyz(1:3,k) =  xyz(1:3,n4(iele,k))
   end do
   call caldiverr(elm_xyz,v3,v3threshold,ediv,denominator,v3abs,iflag)
   call calroterr(elm_xyz,v3,v3threshold,erot,denominator,v3abs,jflag)
!   if ( iflag .eq. 1)   diverr(i) = max(ediv,diverr(i))
   diverr(i) = max(ediv,diverr(i))
   roterr(i) = max(erot,roterr(i))
  end do
 end do

!#[2]## output diverr xy_tri for it
 if ( it .eq. 1) then
 outputfile=outputfolder(1:len_trim(outputfolder))//"diverr_xy.dat"
 open(52,file=outputfile)
  write(52,'(i6)') ntris
  do i=1,ntris
   write(52,'(2g15.7)') xy_tri(1:2,i)
  end do
 close(52)
 end if

!#[3]## out diverr
 write(num,'(i6.6)') it
 outputfile=outputfolder(1:len_trim(outputfolder))//"diverr"//num//".dat"
 open(52,file=outputfile)
  write(52,'(i6)') ntris
  write(52,'(g15.7)') (diverr(i),i=1,ntris)
 close(52)

!#[4]## out roterr
 outputfile=outputfolder(1:len_trim(outputfolder))//"roterr"//num//".dat"
 open(52,file=outputfile)
  write(52,'(i6)') ntris
  write(52,'(g15.7)') (roterr(i),i=1,ntris)
 close(52)

!#[3]##
  write(*,*) "### OUTDIVROTUDISTRIBUTION END!! ###"

return
end
!#################################
!# coded on 2017.11.02
!# iflag = -2 : denominator is too small
!# iflag = -1 : v3abs < v3threshold
!# iflag =  1 : successful
subroutine caldiverr(elm_xyz,v3,v3threshold,ediv,denominator,v3abs,iflag)
implicit none
real(8) ,  intent(in)   :: elm_xyz(3,4)
real(8) ,  intent(in)   :: v3(3,4)
real(8),   intent(in)   :: v3threshold
real(8),   intent(out)  :: denominator
integer(4),intent(out)  :: iflag
real(8),   intent(out)  :: ediv
real(8),   intent(out)  :: v3abs
integer(4)              :: l1,l2,k
real(8)                 :: dx3v3(3,3)
real(8)                 :: gn(3,4)
real(8)                 :: v

iflag = 0
!#[1]## prepare dx3v3
 dx3v3 = 0.d0
 call gradnodebasisfun(elm_xyz,gn,v) ! gn(3,4) [1/km]
 do l1=1,3  ! dx,dy,dz
  do l2=1,3 ! vx,vy,vz
   do k=1,4 ! 4 nodes
    dx3v3(l1,l2) = dx3v3(l1,l2) + gn(l1,k) * v3(l2,k) ! dx3v3 ! [mm/sec/km]=[10^-6 m/sec/m]
   end do
  end do
 end do

!#[2]## prepare v3abs
 v3abs = sqrt((sum(v3(1,:))/4.)**2. + (sum(v3(2,:))/4.)**2. + (sum(v3(3,:))/4.)**2.)

!#[3]## denominator
!#  denominator = |dxvx|^2 + |dyvy|^2 + |dzvz|^2
 denominator = dx3v3(1,1)**2.d0 + dx3v3(2,2)**2.d0 + dx3v3(3,3)**2.d0

!#[4]## calculation
 if ( v3abs .lt. v3threshold       ) then
  iflag = -1
  ediv  = -1.0
  return
 end if

 if ( denominator .lt. 1.d-10 ) then
  iflag = -2 ! denominator is too small
  ediv  = -1.0
  return
 else
  iflag = 1  ! calculation is finished in standard manner
  ediv  = (( dx3v3(1,1) + dx3v3(2,2) + dx3v3(3,3) )**2.d0)/denominator
  return
 end if

end
!#################################
!# coded on 2017.11.02
!# iflag = -2 : denominator is too small
!# iflag = -1 : v3abs < v3threshold
!# iflag =  1 : successful
subroutine calroterr(elm_xyz,v3,v3threshold,erot,denominator,v3abs,iflag)
implicit none
real(8) ,  intent(in)   :: elm_xyz(3,4)
real(8) ,  intent(in)   :: v3(3,4)
real(8),   intent(in)   :: v3threshold
real(8),   intent(out)  :: denominator
integer(4),intent(out)  :: iflag
real(8),   intent(out)  :: erot
real(8),   intent(out)  :: v3abs
integer(4)              :: l1,l2,k
real(8)                 :: dx3v3(3,3)
real(8)                 :: gn(3,4)
real(8)                 :: v

iflag = 0
!#[1]## prepare dx3v3
 dx3v3 = 0.d0
 call gradnodebasisfun(elm_xyz,gn,v) ! gn(3,4) [1/km]
 do l1=1,3  ! dx,dy,dz
  do l2=1,3 ! vx,vy,vz
   do k=1,4 ! 4 nodes
    dx3v3(l1,l2) = dx3v3(l1,l2) + gn(l1,k) * v3(l2,k) ! dx3v3 ! [mm/sec/km]=[10^-6 m/sec/m]
   end do
  end do
 end do

!#[2]## prepare v3abs
 v3abs = sqrt((sum(v3(1,:))/4.)**2. + (sum(v3(2,:))/4.)**2. + (sum(v3(3,:))/4.)**2.)

!#[3]## denominator
!#  denominator = |dxvx|^2 + |dyvy|^2 + |dzvz|^2
  denominator = dx3v3(2,3)**2.d0 + dx3v3(3,2)**2.d0 &
	   &    + dx3v3(3,1)**2.d0 + dx3v3(1,3)**2.d0 &
	   &    + dx3v3(1,2)**2.d0 + dx3v3(2,1)**2.d0

!#[4]## calculation
 if ( v3abs .lt. v3threshold       ) then
  iflag = -1
  erot  = -1.0
  return
 end if

 if ( denominator .lt. 1.d-10 ) then
  iflag = -2 ! denominator is too small
  erot  = -1.0
  return
 else
  iflag = 1  ! calculation is finished in standard manner
  erot  = ( (dx3v3(2,3) - dx3v3(3,2))**2.d0 + (dx3v3(3,1) - dx3v3(1,3))**2.d0 &
   &  +   (dx3v3(1,2) - dx3v3(2,1))**2.d0 )/denominator
  return
 end if

end

!#################################
!# coded on 2017.10.25
subroutine caldivroterror(em_mesh,h_ocean,e_divrot)
implicit none
type(mesh),      intent(in)  :: em_mesh
type(ocean_data),intent(in)  :: h_ocean
type(error),     intent(out) :: e_divrot
real(8)                      :: gn(3,4),v
integer(4)                   :: ntet, node
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:,:) :: n4flag
real(8),   allocatable,dimension(:,:) :: xyz,vxyz
real(8)                      :: elm_xyz(3,4)
real(8)                      :: erot,ediv, denominator
real(8)                      :: erot_total,erot_max,x3_rot(3)
real(8)                      :: ediv_total,ediv_max,x3_div(3)
real(8)                      :: dx3v3(3,3),v3(3,4),v3abs
real(8)                      :: dx3v3_rot(3,3),dx3v3_div(3,3)
real(8)                      :: gn_rot(3,4),gn_div(3,4)
real(8)                      :: v3_rot(3,4),v3_div(3,4)
integer(4)                   :: i,j,l1,l2,k,ii,jj
real(8)                      :: vv(3),a4(4),n3(3) ! for face

!#[1]# set
ntet   = em_mesh%ntet
node   = em_mesh%node
allocate(xyz(3,node),n4(ntet,4),n4flag(ntet,2))
xyz    = em_mesh%xyz
n4     = em_mesh%n4
n4flag = em_mesh%n4flag
allocate(vxyz(3,node))
vxyz   = h_ocean%vxyz    ! [mm/sec]

!#[2]# cal e_rot
erot_total = 0.d0 ; erot_max = 0.d0
ediv_total = 0.d0 ; ediv_max = 0.d0

!============================================================= ocean element loop
do i=1,ntet
 if ( n4flag(i,2) .ne. 2 ) cycle ! when element is not in cean
 !# set v3(3,4) and elm_xyz(3,4)
 do j=1,4
  v3(     1:3,j) = vxyz(1:3,n4(i,j))
  elm_xyz(1:3,j) =  xyz(1:3,n4(i,j))
 end do
 v3abs = sqrt((sum(v3(1,:))/4.)**2. + (sum(v3(2,:))/4.)**2. + (sum(v3(3,:))/4.)**2.)
 if ( v3abs .lt. 10.d0 ) cycle ! 10 mm/sec = 1cm/sec
 !# set gn(3,4)
 call gradnodebasisfun(elm_xyz,gn,v) ! gn(3,4) [1/km]

 dx3v3 = 0.d0
 do l1=1,3  ! dx,dy,dz
  do l2=1,3 ! vx,vy,vz
   do k=1,4 ! 4 nodes
    dx3v3(l1,l2) = dx3v3(l1,l2) + gn(l1,k) * v3(l2,k) ! dx3v3 ! [mm/sec/km]=[10^-6 m/sec/m]
   end do
  end do
 end do

 ! denominator = |dxvx|^2 + |dyvy|^2 + |dzvz|^2
 denominator = dx3v3(1,1)**2.d0 + dx3v3(2,2)**2.d0 + dx3v3(3,3)**2.d0
 if ( denominator .lt. 1.e-10 ) cycle
 ediv  = (( dx3v3(1,1) + dx3v3(2,2) + dx3v3(3,3) )**2.d0)/denominator

 !# another method using face integral for div U
  if ( .false. ) then
   ediv = 0
   call area4(elm_xyz,a4)
   do ii=1,4 ! face loop
    vv(1:3) = (v3(:,lmn(ii,1)) + v3(:,lmn(ii,2)) + v3(:,lmn(ii,3)))/3.d0
    n3(1:3) = -gn(1:3,ii)/sqrt(gn(1,ii)**2. + gn(2,ii)**2. + gn(3,ii)**2.)
    ediv  = ediv + inner(n3,vv)*a4(ii)/v ! [ mm/sec *(km)^2/[km]^3]=[10^-6 m/sec/m]
   end do
   ediv = (ediv**2.d0)/denominator
  end if

 denominator = dx3v3(2,3)**2.d0 + dx3v3(3,2)**2.d0 &
	   &   + dx3v3(3,1)**2.d0 + dx3v3(1,3)**2.d0 &
	   &   + dx3v3(1,2)**2.d0 + dx3v3(2,1)**2.d0
 if ( denominator .lt. 1.e-10 ) cycle
 erot = ( (dx3v3(2,3) - dx3v3(3,2))**2.d0 + (dx3v3(3,1) - dx3v3(1,3))**2.d0 &
   &  +   (dx3v3(1,2) - dx3v3(2,1))**2.d0 )/denominator

 !sum up
 erot_total = erot_total + erot
 ediv_total = ediv_total + ediv
 if ( erot_max .lt. erot ) then
  erot_max  = erot
  x3_rot(1) = sum(elm_xyz(1,:))/4.d0
  x3_rot(2) = sum(elm_xyz(2,:))/4.d0
  x3_rot(3) = sum(elm_xyz(3,:))/4.d0
  dx3v3_rot = dx3v3
  gn_rot    = gn
  v3_rot    = v3
 end if
 if ( ediv_max .lt. ediv ) then
  ediv_max  = ediv
  x3_div(1) = sum(elm_xyz(1,:))/4.d0
  x3_div(2) = sum(elm_xyz(2,:))/4.d0
  x3_div(3) = sum(elm_xyz(3,:))/4.d0
  dx3v3_div = dx3v3
  gn_div    = gn
  v3_div    = v3
 end if
end do

!#[3]## set
!# rot
e_divrot%x3_rot     = x3_rot
e_divrot%erot_total = erot_total
e_divrot%erot_max   = erot_max
e_divrot%dx3v3_rot  = dx3v3_rot
e_divrot%gn_rot     = gn_rot
e_divrot%v3_rot     = v3_rot
!# div
e_divrot%x3_div     = x3_div
e_divrot%ediv_total = ediv_total
e_divrot%ediv_max   = ediv_max
e_divrot%dx3v3_div  = dx3v3_div
e_divrot%gn_div     = gn_div
e_divrot%v3_div     = v3_div

return
end

!######################################################
!# 2017.10.26
!# copied from fluidity-4.1.10/femtools/Funits.F90
  function free_unit()
    !!< Find a free unit number. Start from unit 10 in order to ensure that
    !!< we skip any preconnected units which may not be correctly identified
    !!< on some compilers.
    integer :: free_unit
    
    logical :: connected

    do free_unit=10, 99

       inquire(unit=free_unit, opened=connected)

       if (.not.connected) return

    end do
    
    write(*,*) "No free unit numbers avalable"
    stop

  end function

!###################################################################
!# coded on 2017.10.26
subroutine OPENERRFILE(g_param,e_divrot)
implicit none
type(error),        intent(inout) :: e_divrot
type(param_forward),intent(in)    :: g_param
character(50)                     :: outputfolder, outputfile
integer(4)                        :: idev_div,idev_rot

!#[0]##
outputfolder = g_param%outvxyzfolder ! 2018.11.14

!#[1]## open files
! rot error
outputfile = outputfolder(1:len_trim(outputfolder))//"error_rot.dat"
idev_rot   = free_unit()
open(idev_rot,file=outputfile)

! dev error
outputfile = outputfolder(1:len_trim(outputfolder))//"error_div.dat"
idev_div   = free_unit()
open(idev_div,file=outputfile)

!#[2]## set out device number
e_divrot%idev_rot = idev_rot
e_divrot%idev_div = idev_div

return
end

!###################################################################
!# 2017.10.26
subroutine CLOSEERRFILE(e_divrot)
implicit none
type(error),intent(in) :: e_divrot
integer(4)             :: idev_div,idev_rot

!#[1]##
idev_div = e_divrot%idev_div
idev_rot = e_divrot%idev_rot

!#[2]##
close(idev_div)
close(idev_rot)

return
end
!###################################################################
!# 2017.10.26
subroutine OUTEDIVROTERROR(e_divrot,it,dt)
implicit none
integer(4), intent(in) :: it       ! time step
real(8),    intent(in) :: dt       ! dt [sec]
type(error),intent(in) :: e_divrot
integer(4)             :: idev_rot,idev_div
real(8)                :: x3_rot(3), erot_max,erot_total
real(8)                :: x3_div(3), ediv_max,ediv_total
real(8)                :: time
integer(4)             :: i,j
real(8)                :: dx3v3_rot(3,3),dx3v3_div(3,3)
real(8)                :: gn_rot(3,4),gn_div(3,4)
real(8)                :: v3_rot(3,4),v3_div(3,4)

time = it * dt /60.d0 ! 2017.10.26

!#[1]## output rot error
x3_rot     = e_divrot%x3_rot
erot_total = e_divrot%erot_total
erot_max   = e_divrot%erot_max
idev_rot   = e_divrot%idev_rot
dx3v3_rot  = e_divrot%dx3v3_rot
gn_rot     = e_divrot%gn_rot
v3_rot     = e_divrot%v3_rot
write(idev_rot,'(i7,6g15.7)') it, time, erot_total,erot_max, x3_rot(1:3)
write(idev_rot,*) "dx3v3"
write(idev_rot,'(3g15.7)')  ((dx3v3_rot(i,j),j=1,3),i=1,3)
write(idev_rot,*) "gn"
write(idev_rot,'(4g15.7)')  ((gn_rot(i,j),j=1,4),i=1,3)
write(idev_rot,*) "v3"
write(idev_rot,'(4g15.7)')  ((v3_rot(i,j),j=1,4),i=1,3)

!#[2]## output div error
x3_div     = e_divrot%x3_div
ediv_total = e_divrot%ediv_total
ediv_max   = e_divrot%ediv_max
idev_div   = e_divrot%idev_div
dx3v3_div  = e_divrot%dx3v3_div
gn_div     = e_divrot%gn_div
v3_div     = e_divrot%v3_div
write(idev_div,'(i7,6g15.7)') it, time, ediv_total,ediv_max, x3_div(1:3)
write(idev_div,*) "dx3v3"
write(idev_div,'(3g15.7)')  ((dx3v3_div(i,j),j=1,3),i=1,3)
write(idev_div,*) "gn"
write(idev_div,'(4g15.7)')  ((gn_div(i,j),j=1,4),i=1,3)
write(idev_div,*) "v3"
write(idev_div,'(4g15.7)')  ((v3_div(i,j),j=1,4),i=1,3)

return
end

end module divroterror
