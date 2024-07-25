!# coded by Takuto MINAMI on Aug 31, 2018
module global_ocean_cond_map
implicit none

type ocean_cond
 integer(4)                              :: nlon,nlat,ndepth
 real(8),   allocatable,dimension(:)     :: lon         ! standard lon
 real(8),   allocatable,dimension(:)     :: lat         ! standard lat
 real(8),   allocatable,dimension(:)     :: depth       ![m] standard depth
 real(8),   allocatable,dimension(:,:)   :: lonbnds     ![deg] lon bin
 real(8),   allocatable,dimension(:,:)   :: latbnds     ![deg] lat bin
 real(8),   allocatable,dimension(:,:)   :: depthbnds   ![m]   depth bin (downward
 real(8),   allocatable,dimension(:,:,:) :: val    ! [S/m] global cond grid data
 logical,   allocatable,dimension(:,:,:) :: exist  ! whether cond value exists or not
end type
!
contains


!#####################################################   read global cond
 subroutine read_global_cond(g_cond,woacsvfile)
 implicit none
 type(ocean_cond),intent(out)            :: g_cond
 character(50),   intent(in)             :: woacsvfile
 ! internal variables
 real(8), allocatable,dimension(:)       :: array       ! work array
 real(8),         dimension(2)           :: latlon
 integer(4)                              :: i,iunit,junit,ilen,ii,ilon,ilat
 character(623)                          :: a
 character(3)                            :: num
 real(8),   allocatable,dimension(:)     :: lon        ! standard lon
 real(8),   allocatable,dimension(:)     :: lat        ! standard lat
 real(8),   allocatable,dimension(:)     :: depth      ![m] standard depth
 real(8),   allocatable,dimension(:,:)   :: lonbnds    ![deg] lon bounds for bin
 real(8),   allocatable,dimension(:,:)   :: latbnds    ![deg] lat bounds for bin
 real(8),   allocatable,dimension(:,:)   :: depthbnds  ![m] depth bin(down positive)
 real,      allocatable,dimension(:,:,:) :: val        ! [S/m] global cond grid data
 logical,   allocatable,dimension(:,:,:) :: exist      ! whether cond value exists or not

!# size for woa conductivity 2018.08.31
 integer(4),parameter                   :: nlon    = 360
 integer(4),parameter                   :: nlat    = 180
 integer(4),parameter                   :: ndepth  = 102
 real(8),   parameter                   :: fillval = -9999. ! for abscent bin

!#[0]## set and initialize
  allocate( lon(nlon),lat(nlat),depth(ndepth))
  allocate( lonbnds(2,nlon),latbnds(2,nlat),depthbnds(2,ndepth))
  allocate( val(nlon,nlat,ndepth),exist(nlon,nlat,ndepth))
  allocate( array(ndepth))
  val    = fillval
  exist  = .false.

!#[1]## prepare grids for woa file
!# mean lon, west/east bounds
  do i=1,nlon
   lon(i)       = -180.d0 + 0.5d0 + 1.d0*(i-1) ! [deg]
   lonbnds(1,i) = lon(i) - 0.5d0 ! [deg] west bnds
   lonbnds(2,i) = lon(i) + 0.5d0 ! [deg] east bnds
  end do
!# mean lat, lat west/east bounds
  do i=1,nlat
   lat(i)       = -90.d0 + 0.5d0 + 1.d0*(i-1) ! [deg]
   latbnds(1,i) = lat(i) - 0.5d0 ! [deg] west bnds
   latbnds(2,i) = lat(i) + 0.5d0 ! [deg] east bnds
  end do
!# standard depth
   depth(1:ndepth) =&
 &(/0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,350,375,&
 &400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,&
 &1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,&
 &3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,&
 &5100,5200,5300,5400,5500/) ! [m]
 !# depth bin
  depthbnds(1:2,1)=(/0.d0,2.5d0/)
  do i=2,ndepth - 1
   depthbnds(1,i) = depthbnds(2,i-1)              ! upper depth bound
   depthbnds(2,i) = (depth(i) + depth(i+1))/2.d0  ! lower depth bound
  end do
  depthbnds(1,ndepth) = depthbnds(2,ndepth-1)
  depthbnds(2,ndepth) = depth(ndepth)

!#[2]## read file
  open(newunit=iunit,file=woacsvfile,status='old')
!  open(newunit=junit,file='test.out',status='replace') ! for check
  read(iunit,*) ! header 1
  read(iunit,*) ! header for std depth
  i=0
  do while ( i == 0 )
   array = -1.; a=""
   read(iunit,'(a)',end=99) a ! read as a string
   ilen=len_trim(a)           ! length of the string
   ii=(ilen-11)/6             ! number of cond data for this lat-lon bin
   backspace(iunit)
   read(iunit,*) latlon(1:2),array(1:ii) ! read as real values

   !# search lon lat bin
   call search(nlon,nlat,lonbnds,latbnds,latlon(2),latlon(1),ilon,ilat)
   val(ilon,ilat,1:ii)   = array(1:ii)
   exist(ilon,ilat,1:ii) = .true.

   !# check by output
   !  write(junit,'(a)') a
   write(num,'(i3.3)') 2+ii ! latlon + cond data
!   write(junit,'('//num(1:len_trim(num))//'f10.3)') latlon(1:2), array(1:ii)! for check

  end do
  99 continue
  close(iunit)
!  close(junit) ! for check
  write(*,*) woacsvfile(1:len_trim(woacsvfile))," read end!!"

 !#[3]## output
  g_cond%nlon      = nlon
  g_cond%nlat      = nlat
  g_cond%ndepth    = ndepth
  g_cond%lon       = lon
  g_cond%lat       = lat
  g_cond%depth     = depth
  g_cond%lonbnds   = lonbnds
  g_cond%latbnds   = latbnds
  g_cond%depthbnds = depthbnds
  g_cond%val       = val
  g_cond%exist     = exist

return
end
!#################################################################
subroutine search(nlon,nlat,lonbnds,latbnds,lon,lat,ilon,ilat)
implicit none
real(8),                  intent(in)  :: lon,lat
integer(4),               intent(in)  :: nlon,nlat
real(8),dimension(2,nlon),intent(in)  :: lonbnds
real(8),dimension(2,nlat),intent(in)  :: latbnds
integer(4),               intent(out) :: ilon,ilat
integer(4) :: i,j

do i=1,nlat
 if (latbnds(1,i) .lt. lat .and. lat .lt. latbnds(2,i) ) then
 ilat = i
 goto 100
 end if
end do
write(*,*) "GEGEGE lon not found!"
stop
100 continue

do i=1,nlon
if (lonbnds(1,i) .lt. lon .and. lon .lt. lonbnds(2,i) ) then
 ilon = i
 goto 101
 end if
end do
write(*,*) "GEGEGE lat not found!"
stop
101 continue

return
end

end module

