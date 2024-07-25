! coded on 2018.11.13
module caloceancond
use global_ocean_cond_map ! m_global_ocean_cond.f90
use mesh_type             ! m_mesh_type.f90
use param                 ! m_param.f90
use spherical             ! m_spherical.f90
implicit none

contains
!################################################################### setoceancond
!# coded on 2018.11.13
subroutine setoceancond(g_mesh,g_param)
implicit none
type(mesh),            intent(inout)  :: g_mesh  ! m_mesh_type.f90
type(param_forward),   intent(in)     :: g_param ! m_param.f90
type(ocean_cond)                      :: g_cond  ! m_global_ocean_cond.f90
integer(4)                            :: ntet
integer(4)                            :: node
real(8),   allocatable,dimension(:,:) :: xyz
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:)   :: igroup
real(8),   allocatable,dimension(:)   :: lon,lat,depth,cond
real(8),   allocatable,dimension(:)   :: cmodel
character(50)                         :: woafile
integer(4)                            :: ntet_ocean, nshift, i,j,k
real(8)                               :: lonorigin,latorigin,lonlatalt(3)
real(8),   allocatable,dimension(:,:) :: xc            ! 2018.11.15
integer(4)                            :: iflag_oceancond
integer(4),allocatable,dimension(:)   :: ocean2ntetptr ! 2018.11.14
real(8)                               :: oceancond_default ! 2018.11.15

!# set
 ntet      = g_mesh%ntet
 node      = g_mesh%node
 allocate( igroup(ntet),n4(ntet,4),xyz(3,node) )
 xyz       = g_mesh%xyz
 n4        = g_mesh%n4
 igroup    = g_mesh%n4flag(:,2)
 woafile   = g_param%woafile
 lonorigin = g_param%lonorigin
 latorigin = g_param%latorigin
 iflag_oceancond = g_param%iflag_oceancond
 oceancond_default = g_param%cond(2)  ! 2018.11.15

!# gen default cmodel
 allocate(cmodel(ntet))
 do i=1,ntet
  if ( igroup(i) .eq. 1 ) cmodel(i) = g_param%cond(1) ! air conductivity
  if ( igroup(i) .eq. 2 ) cmodel(i) = g_param%cond(2) ! default ocean cond [S/m]
  if ( igroup(i) .eq. 3 ) cmodel(i) = g_param%cond(3) ! solid earth conductivity
 end do

!# iflag_oceancond = 2 ( woa file is given)
if (iflag_oceancond .eq. 2 ) then ! when woa file is given

 !# read oceancond
  call read_global_cond(g_cond,woafile)

 !# set ocean2ntetptr
  ntet_ocean = 0
  do i=1,ntet
   if ( igroup(i) .eq. 2 ) ntet_ocean = ntet_ocean + 1
  end do
  allocate(ocean2ntetptr(ntet_ocean))
  j=0
  do i=1,ntet
   if ( igroup(i) .eq. 2 ) then
    j=j+1
    ocean2ntetptr(j) = i
   end if
  end do

 !# prepare lon, lat, depth
  allocate( lon(ntet_ocean),lat(ntet_ocean),depth(ntet_ocean),xc(3,ntet_ocean))
   xc = 0.d0
  do i=1,ntet_ocean
   do k=1,4
    xc(1:3,i)=xc(1:3,i) + xyz(:,n4(ocean2ntetptr(i),k))/4.d0
   end do
   call xyz2lonlatalt(xc(:,i),lonorigin,latorigin,lonlatalt)
   lon(i)   =   lonlatalt(1)
   lat(i)   =   lonlatalt(2)
   depth(i) = - lonlatalt(3)*1000.
  end do

 !# cal ocean cond
  allocate( cond(ntet_ocean))
  call calpointoceancond(ntet_ocean,lon,lat,depth,cond,woafile)
  open(1,file="oceancond.dat")
  do i=1,ntet_ocean
   if ( cond(i) .le. 0. ) cond(i) = oceancond_default ! 2018.11.15 avoid -9999.0, for empty bin
   cmodel(ocean2ntetptr(i)) = cond(i)
   write(1,'(i10,g15.7,a,3f12.5)') ocean2ntetptr(i),cond(i),"xyz",xc(1:3,i)
  end do
  close(1)
 end if

!# generate output
  g_mesh%cmodel = cmodel

return
end

!####################################################################  caloceancond
!# coded on Aug 31, 2018
 subroutine calpointoceancond(npoint,lon,lat,depth,cond,woafile)
 use global_ocean_cond_map
 implicit none
 integer(4),            intent(in)       :: npoint
 real(8),               intent(in)       :: lon(npoint),lat(npoint),depth(npoint)
 real(8),               intent(out)      :: cond(npoint)
 character(50),         intent(in)       :: woafile
 type(ocean_cond)                        :: g_cond
 integer(4)                              :: i,j,k,ii,ilon,jlat,kdepth
 integer(4)                              :: nlon,nlat,ndepth
 real(8),   allocatable,dimension(:,:)   :: lonbnds  ![deg] lon bin
 real(8),   allocatable,dimension(:,:)   :: latbnds  ![deg] lat bin
 real(8),   allocatable,dimension(:,:)   :: depthbnds! [m] depth bin(down positive)
 real(8),   allocatable,dimension(:,:,:) :: val      ! [S/m] global cond grid data
 real(8)                                 :: dd       ! 2018.11.14

!#[1]# get global ocean cod
 call read_global_cond(g_cond,woafile)


!#[2]# set ##
 nlon      = g_cond%nlon
 nlat      = g_cond%nlat
 ndepth    = g_cond%ndepth
 lonbnds   = g_cond%lonbnds
 latbnds   = g_cond%latbnds
 depthbnds = g_cond%depthbnds

!#[3]## cal cond(1:npoint)
  do ii=1,npoint

   !# depth adjust 2018.11.14
   dd = depth(ii)
   if ( depth(ii) .gt. depthbnds(2,ndepth) ) then
    dd = depthbnds(2,ndepth)
   end if

   do i=1,nlon
    if ( lonbnds(1,i) .le. lon(ii) .and. lon(ii) .le. lonbnds(2,i)) then
    do j=1,nlat
     if ( latbnds(1,j) .le. lat(ii) .and. lat(ii) .le. latbnds(2,j)) then
     do k=1,ndepth
     if ( depthbnds(1,k) .le. dd .and. dd .le. depthbnds(2,k)) then
      ilon   = i
      jlat   = j
      kdepth = k
      goto 100
     end if
     end do
     end if
    end do
    end if
   end do
   write(*,*) "GEGEGE not found lon lat depth=",lon(ii),lat(ii),dd
   stop
100 continue
    cond(ii) = g_cond%val(ilon,jlat,kdepth)
  end do

return
end


end module


