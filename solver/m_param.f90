! Coded on T.M. on May 16, 2016
module param
implicit none

type param_ana ! parameters for analytical velocity
  real(8)      :: depth       ![km]  ocean depth for linear long wave
  real(8)      :: lambda      ![km]  wavelength of tsunami
  real(8)      :: eta         ![m]   tsunami height amplitude
  real(8)      :: tlarge      ![sec] time to full wave height
  real(8)      :: xwavelength ![1  ] how many wavelength in x direction
  real(8)      :: yextent     ![km]  horizontal extent of wave in y direction
end type param_ana

type param_forward
 integer(4)    :: icalerrflag !1:err,2:forward calculation, 3:both 2017.10.25
 character(70) :: em3dmsh     ! 3D EM mesh file
 character(70) :: linefile    ! line info of 3D EM mesh
 character(70) :: meshctlfile ! mesh control file
 character(70) :: oceanfile  ! ocean.msh for nodes

 !#[2]## geomag info
 integer(4)    :: iflag_geomag ! 1 for file, 2 for static vector
 !#[2-0]## iflag_geomag = 0  ! IGRF12 is used 2018.11.14
 real(8)       :: year_decimal ! decimal year for IGRF12, grid is same as
 !#[2-1]## iflag_geomag = 1
 character(70) :: xgrdfile_f ! eastgrid (used for iflag_geomag = 0 or 1)
 character(70) :: ygrdfile_f ! northgrid
 character(70) :: geomagfile ! geomag file (grid should be the same as M, N)
 !#[2-2]## iflag_geomag = 2
 real(8)       :: geomagvector(1:3) ! east,north,upward [nT]
 integer(4)    :: iflag_vF   ! 1:vh*Fz, 2:vz*Fh, 3: full (vh*Fz + vz*Fh) 2017.07.03

 !
 character(70) :: outbxyzfolder ! 2018.11.14
 character(70) :: outexyzfolder ! 2018.11.14
 character(70) :: outvxyzfolder ! 2018.11.14

 !#[3]## water velocity info
 integer(4)    :: iflag_velocity ! 1 for comcot, 2 for analytic
 !#[3-1]## COMCOT (iflag_velocity = 1 )
 character(70) :: xgrdfile_v ! eastgrid
 character(70) :: ygrdfile_v ! northgrid
 character(70) :: hgrdfile   ! ocean depth file
 character(70) :: mpre       ! header for M file from comcot
 character(70) :: npre       ! header for N file from comcot
 real(8)       :: dt_comcot  ! time interval of comcot result
 !#[3-2]## analytical velocity (iflag_velocity = 2)
 type(param_ana) :: h_anapara

 !#[4]## time parameters
 real(8)       :: tstart ! [sec]
 real(8)       :: dt     ! [sec]
 real(8)       :: tmax   ! [sec]

 !#[5]## observatory info
 integer(4)    :: nobs  ! # of observatories
 integer(4)    :: lonlatflag !# 1: lonlat [deg], 2: xy [km]
 real(8)       :: lonorigin ! [deg] get from g_meshpara on 2016.11.20
 real(8)       :: latorigin ! [deg] get from g_meshpara on 2016.11.20
 real(8),allocatable,dimension(:,:)     :: lonlataltobs! lon [deg], lat [deg], alt [km]
 real(8),allocatable,dimension(:,:)     :: xyzobs      ! x [km], y[km],z [km] (upward)
 character(70),allocatable,dimension(:) :: obsname

 !#[6]## output option
 integer(4)    :: ixyflag
 ! ixyflag = 0 no xy values are calculated
 ! ixyflag = 1 nx * ny grid surface view
 ! ixyflag = 2 xy values for all the surface triangles
 integer(4)    :: ixydepth
 ! ixydepth = 1 means same depth values are calculated
 ! ixydepth = 2 means values 1m beneath the seafloor/ground surface
 !#[6-0]## when ixyflag = 0
 !# no further parameter input
 !#[6-1]## when ixyflag = 1
 real(8)       :: depth_xy
 integer(4)    :: nx,ny     ! grid for same depth value
 real(8)       :: tinterval_xy ! [sec] time interval for xy plain output 2018.11.14 for ixyflag=1,2
 !#[6-2]## when ixyflag = 2, extrude mesh data are required
 character(70) :: pokimsh   ! polygon_ki mesh file
 character(70) :: ki23dfile! node pointer from pokimsh to 2D ocean surface id

 !#[7]## conductivity information
 ! normally, cond(1) for air, cond(2) for ocean, cond(3) for crust
 integer(4) :: ncond
 real(8),allocatable,dimension(:)  :: cond

 !#[8]##
 integer(4) :: iflagspherical ! 0 for normal, 1 for spherical coordinate
 real(8)    :: xyzminmax(6)

!#[9]## Inclusion of global grid conductivity  2018.08.30
 integer(4)    :: iflag_oceancond ! 1: fixed (or not given), 2: give global cond data
 character(70) :: woafile      ! name global conductivity data from WOA
 end type


contains
!################################################## subroutine readparam
subroutine readparam(c_param)
implicit none
type(param_forward),intent(out) :: c_param
integer(4) :: i,j,input=5

!open(input,file="tohoku.ctl")
read(input,*) ! header
!#[0]## read mesh_param file
write(*,*) "Please input 1 for vxyz error, 2:forward calculation, 3:both" ! 2017.10.25
read(input,12) c_param%icalerrflag    ! 2017.10.25
write(*,*) "icalerrflag =",c_param%icalerrflag ! 2017.10.26
write(*,*) "Please enter em3d meshfile"
read(input,10) c_param%em3dmsh ! line 2
write(*,*) "Please enter lineinfo file"
read(input,10) c_param%linefile
write(*,*) "Please enter mesh control file"
read(input,10) c_param%meshctlfile
write(*,*) "Please enter ocean.msh file path"
read(input,10) c_param%oceanfile
write(*,*) "oceanfile =",c_param%oceanfile ! 2017.10.26

!#[1]## time parameters
write(*,*) "Please enter start time, time interval, end time in seconds"
read(input,11) c_param%tstart
read(input,11) c_param%dt
read(input,11) c_param%tmax
write(*,*)   "tstart,dt,dmax=",c_param%tstart,c_param%dt,c_param%tmax

!#[2]## info for seawater velocity
write(*,*) "Please choose 1 from 2 modes for ocean velocity field"
write(*,*) "1: provide by COMCOT file, 2: analytical velocity"
read(input,12) c_param%iflag_velocity
write(*,*) "iflag_velocity=",c_param%iflag_velocity ! 2017.10.26
 !#[2-1]## COMCOT
if ( c_param%iflag_velocity .eq. 1 ) then
 write(*,*) "Please enter x grid file"
 read(input,10) c_param%xgrdfile_v
 write(*,*) "Please enter y grid file"
 read(input,10) c_param%ygrdfile_v
 write(*,*) "Please enter h grid file"
 read(input,10) c_param%hgrdfile
 write(*,*) "Please enter M header file"
 read(input,10) c_param%mpre
 write(*,*) "Please enter N header file"
 read(input,10) c_param%npre
 write(*,*) "Please enter time interval for COMCOT results"
 read(input,11) c_param%dt_comcot
 !#[2-2]## analytical wave
else if (c_param%iflag_velocity .eq. 2 ) then
 write(*,*) "Please input parameters for analytical velocity for LLW"
 write(*,*) "This program deals with only the tsunami propagates eastward."
 write(*,*) "Please enter ocean depth [km]"
  read(input,11) c_param%h_anapara%depth
 write(*,*) "Please enter wavelength [km]"
  read(input,11) c_param%h_anapara%lambda
 write(*,*) "Please enter tsunami amplitude (height) [m]"
  read(input,11) c_param%h_anapara%eta
 write(*,*) "Please enter time to full tsunami heigh (tlarge) [s]"
  read(input,11) c_param%h_anapara%tlarge
 write(*,*) "Please enter how many wavelengths in x direction (real)"
  read(input,11) c_param%h_anapara%xwavelength
 write(*,*) "Please enter how long the northward extent is [km]"
  read(input,11) c_param%h_anapara%yextent
 write(*,*) "depth =",c_param%h_anapara%depth
 write(*,*) "lambda =",c_param%h_anapara%lambda
 write(*,*) "eta =",c_param%h_anapara%eta
 write(*,*) "tlarge =",c_param%h_anapara%tlarge
 write(*,*) "xwavelength =",c_param%h_anapara%xwavelength
 write(*,*) "yextent =",c_param%h_anapara%yextent
! stop
else
 write(*,*) "GEGEGE! iflag_velocity should be 1 or 2"
 stop
end if

!#[2.5]## iflag_vF   2017.10.26
write(*,*) "iflag_vF:1 for vh*Fz, 2 for vz*Fh, 3 for both"
read(*,12) c_param%iflag_vF

!#[3]## geomag field info
write(*,*) "Please input 0, 1, 2 (iflag_geomag) for background geomagnetic field"!2018.11.14
write(*,*) "0: use IGRF12 with COMCOT flow input (not available when using analytic flow)"
write(*,*) "1: provide input geomag file, 2: constant geomag vector"
read(input,12)  c_param%iflag_geomag
write(*,*)  "iflag_geomag=",c_param%iflag_geomag
if (            c_param%iflag_geomag .eq. 0 ) then                                  ! 2018.11.14
 write(*,*) "IGRF will be used with the COMCOT grid info (read as X:lon, Y:lat)"    ! 2018.11.14
 write(*,*) "Please input decimal year for calculation of background magnetic field"! 2018.11.14
 if ( c_param%iflag_velocity .eq. 2 ) then             ! 2018.11.14
  write(*,*) "GEGEGE! analytical velocity and IGRF is not used together" ! 2018.11.14
  stop  ! 2018.11.14
 end if ! 2018.11.14
 read(input,11) c_param%year_decimal                   ! 2018.11.14
 write(*,*) "Please enter x grid file for geomag data" ! 2018.11.14
 read(input,10) c_param%xgrdfile_f                     ! 2018.11.14
 write(*,*) "Please enter y grid file for geomag data" ! 2018.11.14
 read(input,10) c_param%ygrdfile_f                     ! 2018.11.14
else if (       c_param%iflag_geomag .eq. 1 ) then
 write(*,*) "Please enter x grid file for geomag data"
 read(input,10) c_param%xgrdfile_f
 write(*,*) "Please enter y grid file for geomag data"
 read(input,10) c_param%ygrdfile_f
 write(*,*) "Please enter geomag file (north[nT],east[nT],down[nT])"
 read(input,10) c_param%geomagfile
else if (       c_param%iflag_geomag .eq. 2 ) then
 write(*,*) "Please enter vector geomagnetic field (east,north,upward)[nT]"
 read(input,13) c_param%geomagvector(1:3)
else
 write(*,*) "GEGEGE! iflag_geomag should be 1 or 2"
 stop
end if

!#
write(*,*) "Please enter output bxyz, exyz, vxyz folder with the end of /"
read(input,10) c_param%outbxyzfolder ! 2018.11.14
read(input,10) c_param%outexyzfolder ! 2018.11.14
read(input,10) c_param%outvxyzfolder ! 2018.11.14

!#[5]## observatory info
write(*,*) "Please enter # of observatories"
read(input,12) c_param%nobs
write(*,*) "# of observatories is",c_param%nobs
allocate(c_param%lonlataltobs(3,c_param%nobs))
allocate(   c_param%xyzobs(3,c_param%nobs))
allocate(c_param%obsname(c_param%nobs))
read(input,12) c_param%lonlatflag  !# 1: lonlat [deg], 2: xy [km]
! if lonlatflag .ne. 1, lonorigin and latorigin will not be used
!read(input,*) c_param%lonorigin,c_param%latorigin
if ( c_param%lonlatflag .eq. 1) then
 write(*,*) "Input of coordinates are lonlat [deg] and alt [km] (upward)"
 do i=1,     c_param%nobs
  read(input,10) c_param%obsname(i)
  read(input,13) (c_param%lonlataltobs(j,i),j=1,3)
 end do
else if ( c_param%lonlatflag .eq. 2 ) then
 write(*,*) "Input coordinates are xyz [km]"
 do i=1,c_param%nobs
  read(input,10) c_param%obsname(i)
  read(input,13) (c_param%xyzobs(j,i),j=1,3)
!  write(*,*)    (c_param%xyobs(i,j),j=1,2),c_param%zobs(i)
 end do
else
 goto 99
end if

!#[6]## output option
write(*,*) "Input ixyflag, 0 for nothing, 1 for ouput of planview of result"
read(input,12) c_param%ixyflag
write(*,*) "ixyflag =",c_param%ixyflag

!#[6-0]## ixyflag = 0 : no plan view and rigid given obs z coordinate is used
!# ixyflag = 0 no xy values are calculated (no modification for observatory position)
!# ixyflag = 1 nx * ny grid surface view   (modification to observatory mposition)
!# ixyflag = 2 xy values for all the surface triangles (modification for observatory)
!### depth parameter for
!# ixydepth = 1 means same depth values are calculated
!# ixydepth = 2 means values 1m beneath the seafloor/ground surface
if (c_param%ixyflag .eq. 0 ) then
 write(*,*) "parameter read is finished"

!#[6-1]## ixyflag = 1
 else if (c_param%ixyflag .eq. 1) then
  write(*,*) "Input nx and ny of xy plain grid like 'nx ny [Enter]' "
  read(*,14) c_param%nx,c_param%ny     ! grid for same depth value
  write(*,*) "nx=",c_param%nx,"ny=",c_param%ny
  write(*,*) "Input time interval [sec] for xy plain view data"
  read(*,11) c_param%tinterval_xy
 else if (c_param%ixyflag .eq. 2) then
  write(*,*) " Output the files for all the surface triangles"
  write(*,*) "Input time interval [sec] for xy plain view data"
  read(*,11) c_param%tinterval_xy
end if


if (c_param%ixyflag .ge. 1) then
 write(*,*) "parameter for xy depth:"
 write(*,*) "ixydepth =1 for homogeneous depth"
 write(*,*) "ixydepth =2 for 1m beneath the seafloor/ground surface"
 read(*,12)  c_param%ixydepth
 if (c_param%ixydepth .eq. 1) then
  write(*,*) "Input depth [km] for xy plain data"
  read(*,11) c_param%depth_xy
 end if  
end if


!#[6-2]## ixyflag = 1, 2
if ( c_param%ixyflag .ge. 1 ) then
 write(*,*) "Please enter polygonki file name" ! added on Aril 21, 2016
 read(input,10) c_param%pokimsh
 write(*,*) "Please enter ki23dfile name" ! added on Aril 21, 2016
 read(input,10) c_param%ki23dfile
end if

!#[7]## conductivity information
 write(*,*) "Input # of conductivity groups"
 read(input,12) c_param%ncond
 allocate( c_param%cond(c_param%ncond))
 write(*,*) "Input conductivity values for each group [S/m]"
 do i=1,c_param%ncond
  read(input,11) c_param%cond(i)
 end do

!#[8]## global conductivity is used or not 2018.08.30
 c_param%iflag_oceancond = 1 ! 2018.08.30
 write(*,*) "Input 1 for fixed conductivity, 2 for global conductivity data"!2018.08.30
 read(input,12,end=98) c_param%iflag_oceancond                      !2018.08.30
 if ( c_param%iflag_oceancond .eq. 2 ) then                         !2018.08.30
  write(*,*) "Input csv file name for global coean cond from WOA"   !2018.11.13
  read(input,10)       c_param%woafile                              !2018.08.30
 end if                                                             !2018.08.30

98 continue ! 2018.08.30

10 format(20x,a70)
11 format(20x,g15.7)
12 format(20x,i10)
13 format(20x,3g15.7)
14 format(20x,2i10)
15 format(20x,a10) ! 2018.08.30

write(*,*) "## READPARAM END!! ###"
return

99 continue
write(*,*) "GEGEGE!"
write(*,*) "c_param%lonlatflag=",c_param%lonlatflag
stop
end subroutine readparam
!

end module param
