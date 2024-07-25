! coded on May 20, 2016
!
module obs_type
use matrix
implicit none

type obs_info
 integer(4) :: nobs
 character(4) :: name
 real(8),pointer,dimension(:,:)  :: xyz_obs
 real(8),pointer,dimension(:,:)  :: xyzspherical ! 2016.11.20
 real(8),pointer,dimension(:,:)  :: lonlatalt    ! 2016.11.20
 type (real_crs_matrix)          :: coeff(2,3)
 type (real_crs_matrix)          :: coeff_vF
 integer(4),pointer,dimension(:) :: devnum(:)
end type

contains

!################################################## OPENOBSFILE
subroutine OPENOBSFILE(obs_sites,g_param)
use param
implicit none
type(param_forward),intent(in)    :: g_param
type(obs_info),     intent(inout) :: obs_sites
character(70)                     :: filename, obsname
integer(4)                        :: inum, i
character(70)                     :: outbxyz ! 2018.11.14
integer(4)                        :: lenbxyz ! 2018.11.14

allocate(obs_sites%devnum(obs_sites%nobs))
outbxyz = g_param%outbxyzfolder ! 2018.11.14
lenbxyz = len_trim(outbxyz)     ! 2018.11.14

do i=1,obs_sites%nobs
 inum=60+i
 obs_sites%devnum(i)=inum
 obsname = g_param%obsname(i)
 !# 2018.11.14 output folder is modified from bxyz/ to outputfolder
 filename=outbxyz(1:lenbxyz)//obsname(1:len_trim(obsname))//"_bxyz_ts.dat"
 open(inum,file=filename)
end do

write(*,*) "### OPENOBSFILE END!! ###"

return
end

!################################################## CLOSEOBSFILE
subroutine CLOSEOBSFILE(obs_sites)
implicit none
type(obs_info),intent(in) :: obs_sites
integer(4) :: i

do i=1,obs_sites%nobs
 close(obs_sites%devnum(i))
end do

write(*,*) "### CLOSEOBSFILE END!! ###"

return
end

!################################################## EBXYZOUT_TS
subroutine EBXYZOUT_TS(obs_sites,bxyz_sites,it,dt)
implicit none
type(obs_info),intent(in) :: obs_sites
integer(4),intent(in) :: it
real(8),intent(in) :: dt, bxyz_sites(obs_sites%nobs,3)
integer(4) :: i
real(8) :: t,bxyz_ana(3)

t=it*dt
do i=1, obs_sites%nobs
 call calbxyzana(obs_sites%xyz_obs(1:3,i),it,dt,bxyz_ana)
 write(obs_sites%devnum(i),'(7g15.7)') t,bxyz_sites(i,1:3),bxyz_ana(1:3)
end do

write(*,*) "### EBXYZOUT_TS END!! ###"

return
end

!##################################################
subroutine calbxyzana(x3,it,dt,bxyz)
implicit none
real(8),intent(in) :: x3(3),dt
real(8),intent(out) :: bxyz(3)
integer(4),intent(in) :: it
real(8) :: v3(3),Fz=-30000.d0 ! [nT]
real(8),parameter :: h=8000.d0, g=9.8d0
real(8),parameter :: sigma=3.3d0, pi=4.d0*atan(1.d0)
real(8),parameter :: mu=4.*pi*1.d-7
complex(8) :: iunit=(0.d0,1.d0)
complex(8) :: CT,eta,bz
real(8) :: t,k,c,ch
t=(it+0.5d0)*dt

c=sqrt(g*h) ! [m/sec]
ch=2./mu/sigma/h
k=2.*pi/80.d0 ! [rad/km]

eta=cdexp(iunit*(k*(x3(1)-(c*1.d-3)*t))) ! [m]
bz=c/(c+iunit*ch)*eta/h*Fz
bxyz(3)=real(bz)
bxyz(1)=real(iunit*bz)
bxyz(2)=0.d0
return
end

end module obs_type
