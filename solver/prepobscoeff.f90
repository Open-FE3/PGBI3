!####################################################### PREPOBSCOEFF
! Modified by T. MINAMI on 2016.11.13 to search 1m below the seafloor
! copied from ../FEM_node/n_bzfem.f90
! adjusted to edge-FEM code
subroutine PREPOBSCOEFF(em_mesh,l_line,obs_sites,g_param)
use mesh_type
use line_type
use matrix
use obs_type
use fem_edge_util
use param      ! include lonorigin and latorigin on 2016.11.20
use spherical  ! added on 2016.11.19
implicit none

!integer(4),intent(in) :: nobs ! # observatory
type(mesh),           intent(in)  :: em_mesh
type(line_info),      intent(in)  :: l_line
type(obs_info),       intent(out) :: obs_sites
type(param_forward),  intent(in)  :: g_param
type(mesh) :: ki_mesh
integer(4),allocatable,dimension(:,:) :: ki23dptr
character(70) :: ki23dfile
! (1,(1,2,3)) for edge basis fun (x,y,z), (2,(1,2,3)) for face basis fun
type(grid_list_type) :: glist
integer(4) :: nx,ny,nz,nobs,ntet,nodek,ndat,i,j,ii,ixyflag,iflagspherical
real(8)    :: lonorigin,latorigin,lonlatalt(3)
real(8)    :: xyzminmax(6) ! 1 for normal

!#[0]## set input
 ntet           = em_mesh%ntet
 nobs           = g_param%nobs
 ki23dfile      = g_param%ki23dfile
 obs_sites%nobs = nobs
 obs_sites%name = "Site" !up to 4 characters
 allocate(obs_sites%xyz_obs(3,nobs))
 obs_sites%xyz_obs = g_param%xyzobs
 ixyflag = g_param%ixyflag
 !#[OPTION]# spherical calculation
 iflagspherical = g_param%iflagspherical ! added on 2016.11.19
 lonorigin = g_param%lonorigin ! added pn 2016.11.20
 latorigin = g_param%latorigin ! added on 2016.11.20
 xyzminmax      = g_param%xyzminmax

!#[OPTION]## if ixyflag >= 1, modify the z coordinate of observatory
  !#[O-1]## read ki_mesh and set input  added on 2016.11.13
  if ( ixyflag .ge. 1) then ! 2016.11.20
   CALL READMESH_TOTAL(ki_mesh, g_param%pokimsh)
   nodek  = ki_mesh%node
   !#[O-2]## read ki23dptr
   allocate(ki23dptr(2,ki_mesh%node))
   open(1,file=ki23dfile)
    read(1,*) ndat
    if ( ndat .ne. nodek) goto 99
    do i=1,ndat
     read(1,*) ki23dptr(1:2,i)
    end do
   close(1)
   !#[O-3]## modify zcoordinate of obs 1m beneath the seafloor
   CALL PREPZOBS(obs_sites,em_mesh,ki_mesh,ki23dptr,g_param,0) ! get zobs, 0 for nothing
  end if

!#[OPTION]## iflagspherical = 1  ##
  if (iflagspherical .eq. 1) then
    allocate(obs_sites%xyzspherical(3,nobs)) ! added on 2016.11.19
    allocate(obs_sites%lonlatalt(3,nobs)   ) ! added on 2016.11.19
    do i=1,nobs
     call xyz2lonlatalt(obs_sites%xyz_obs(:,i),lonorigin,latorigin,obs_sites%lonlatalt(:,i))
     call lonlatalt2xyzspherical(obs_sites%lonlatalt(:,i),obs_sites%xyzspherical(:,i))
     write(*,*) "obs_sites%lonlatalt(1:3,i)=",obs_sites%lonlatalt(1:3,i)
     write(*,*) "obs_sites%xyzspherical(1:3,i)=",obs_sites%xyzspherical(1:3,i)
    end do
  end if

!#[3]## generate grid and classify all the elements to gridded boxes
 nx=300;ny=300;nz=300 ! please lower if corresponding element was not found
 write(*,*) "ntet=",ntet,"nx,ny,nz=",nx,ny,nz
 CALL allocate_grid_list(nx,ny,nz,ntet,glist)   ! see m_mesh_type.f90
 CALL gengridforlist(xyzminmax,glist)   ! normal
 CALL classifyelement2grd(em_mesh,glist)! m_mesh_type.f90
 write(*,*) "categolize elements end"

!#[4]# allocate coeffobs and calculate
!cal coeffobs(1:2,1:3)
 CALL CALCOEFF_LINE(em_mesh,l_line,glist,obs_sites,iflagspherical)

 write(*,*) "### PREPOBSCOEFF END!! ###"
 return

99 continue
 write(*,*) "GEGEGE! ndat=",ndat,"is not equal to nodek",nodek
 stop
end

!########################################################### PREPZOBS
! Copied from volcano/3D_ana_comp/FEM_edge_bxyz_model/n_ebfem_bxyz.f90
! on 2016.11.13
! calcualte h_mesh%zobs from h_mesh%xyz
subroutine PREPZOBS(obs,em_mesh,h_mesh,ki23dptr,g_param,iflag)
use mesh_type
use obs_type ! 2016.11.20
use matrix   ! 2016.11.20
use triangle ! 2016.11.23 see m_triangle
use param    ! 2016.11.23
use triangle ! 2016.11.23
implicit none
type(mesh),         intent(in)    :: em_mesh,h_mesh
integer(4),         intent(in)    :: ki23dptr(2,h_mesh%node)
type(param_forward),intent(in)    :: g_param
type(obs_info),     intent(inout) :: obs
integer(4),         intent(in)    :: iflag ! 0: nothing, 1 for calculate coeff_vF for IXYHOUT
real(8),            allocatable,dimension(:,:) :: xyzk,xyz ! triangle,tetra nodes
real(8),            allocatable,dimension(:,:) :: xyzobs
integer(4),         allocatable,dimension(:,:) :: n3k
type(grid_list_type) :: glist
integer(4)           :: ntri,nobs
integer(4)           :: i,j,k,n1,n2,n3,m1,m2,m3,nx,ny,iele
real(8),dimension(2) :: x12,x13,x23,v1,v2,v3
real(8) :: a3(3),a,xyzminmax(6)
!#[0]## set
nobs      = obs%nobs
allocate( xyzk(3,h_mesh%node),xyz(3,em_mesh%node) )
allocate( xyzobs(3,nobs),     n3k(h_mesh%ntri,3)  )
xyzk      = h_mesh%xyz    ! triangle mesh
xyz       = em_mesh%xyz   ! tetrahedral mesh
n3k       = h_mesh%n3
ntri      = h_mesh%ntri
xyzobs    = obs%xyz_obs
xyzminmax = g_param%xyzminmax
if (iflag .eq. 1) call allocate_real_crs_with_steady_ncolm(obs%coeff_vF,nobs,3)
write(*,*) "PREPOBSCOEFF start!"

!#[0]## generate element list
nx=300;ny=300
CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_mesh_type.f90
CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see

!#[1] search for the triangle including (x1,y1)
do j=1,nobs

call findtriwithgrid(h_mesh,glist,xyzobs(1:2,j),iele,a3)
! do i=1,ntri
    n1 = n3k(iele,1)
    n2 = n3k(iele,2)
    n3 = n3k(iele,3)
!   x12(1:2) = xyzk(1:2,n2) - xyzk(1:2,n1) ! [km]
!   x13(1:2) = xyzk(1:2,n3) - xyzk(1:2,n1) ! [km]
!   x23(1:2) = xyzk(1:2,n3) - xyzk(1:2,n2) ! [km]
!   v1 = xyzobs(1:2,j) - xyzk(1:2,n1)
!   v2 = xyzobs(1:2,j) - xyzk(1:2,n2)
!   v3 = xyzobs(1:2,j) - xyzk(1:2,n3)
!   a =( x13(1)*x12(2) - x13(2)*x12(1))
!   a2=( x13(1)* v1(2) - x13(2)* v1(1))/a
!   a3=(  v1(1)*x12(2) -  v1(2)*x12(1))/a
!   a1=(  v2(1)*x23(2) -  v2(2)*x23(1))/a
!   if (a1 .ge. 0. .and. a2 .ge. 0. .and. a3 .ge. 0. ) then
    m1 = ki23dptr(2,n1)
    m2 = ki23dptr(2,n2)
    m3 = ki23dptr(2,n3)
    !#[OPTION] for IXYHOUT ########## 2016.11.20
    if ( iflag .eq. 1) then
     obs%coeff_vF%stack(j)=3*j
     obs%coeff_vF%val((j-1)*3+1) = a3(1)
     obs%coeff_vF%val((j-1)*3+2) = a3(2)
     obs%coeff_vF%val((j-1)*3+3) = a3(3)
     obs%coeff_vF%item((j-1)*3+1) = ki23dptr(1,n1) ! surface
     obs%coeff_vF%item((j-1)*3+2) = ki23dptr(1,n2) ! surface
     obs%coeff_vF%item((j-1)*3+3) = ki23dptr(1,n3) ! surface
!     write(*,*) "j=",j,"coeff_vF%val((j-1)*3+1:3)=",(obs%coeff_vF%val((j-1)*3+k),k=1,3)
!     write(*,*) "coeff_vF%item((j-1)*3+1:3)=",(obs%coeff_vF%item((j-1)*3+k),k=1,3)
    end if !############################### 2016.11.20
!    write(*,*) "n1=",n1,"n2=",n2,"n3=",n3
!    write(*,*) "m1=",m1,"m2=",m2,"m3=",m3
    obs%xyz_obs(3,j) = a3(1)*xyz(3,m1)+a3(2)*xyz(3,m2)+a3(3)*xyz(3,m3) - 0.001! 1m below seafloor
    if ( iflag .ne. 1) then ! for real observatories
     write(*,*) "xyobs(1:2,j)=",xyzobs(1:2,j)
     write(*,*) j,"/nobs",xyzobs(3,j),"->",obs%xyz_obs(3,j),"[km]"
    end if
!   write(*,*) "element # =",i
!   write(*,*) "a1,a2,a3=",a1,a2,a3
!   write(*,*) "xyzg(3,n1 - n3)=",xyzg(3,n1),xyzg(3,n2),xyzg(3,n3)
!    goto 100
!   end if
! end do
 100 continue
end do

write(*,*) "### PREPZOBS END!! ###"
return
end

