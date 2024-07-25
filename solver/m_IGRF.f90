!# coded on 2018.11.14
module IGRF
use param
use FROMCOMCOT
implicit none

contains

!#####################################
subroutine CALGRDIGRF(h_grd,g_param)
implicit none
	type(param_forward), intent(in)    :: g_param
	type(grid_data_2D),  intent(inout) :: h_grd
	integer(4) :: isv, nx, ny, i,j, itype
	real(8)    :: date, x, y, z, f, elong, colat, alt
	real(8),allocatable,dimension(:)  :: lon, lat
	character(70) :: outputfolder,outfile ! 2018.11.15


!#[1] set
	write(*,*) "IGRF calculation"
	date = g_param%year_decimal
	write(*,*) "Decimal year for IGRF =", date
	alt=0.d0 ! altitude [km]
	write(*,*) "Altitude [km] =",alt
	nx  = h_grd%xygrid%nx
	ny  = h_grd%xygrid%ny
	allocate(lon(nx),lat(ny))
	lon = h_grd%xygrid%xgrd
	lat = h_grd%xygrid%ygrd
	outputfolder = g_param%outbxyzfolder ! 2018.11.15
	outfile = outputfolder(1:len_trim(outputfolder))//"geomag_IGRF.dat" ! 2018.11.15

!#[3]## calculate FH, Fz and out put
      isv=0    ! main-field values
	itype=1  ! geodetic
      do i=1,ny
	  do j=1,nx
	    elong=lon(j)             ! east-longitude [0-360]
	    colat=90.d0 - lat(i) ! colatitiude [0-180]
          call igrf12syn (isv,date,itype,alt,colat,elong,x,y,z,f) ! x: north, y: east, z: down
          h_grd%dat(1:3,j,i) = (/x,y,z/)
	  end do
	end do
	close(1)

!#[4]##
	open(1,file=outfile)
	 do i=1,ny
	  do j=1,nx
	   write(1,'(5g15.7)') lon(j),lat(i),h_grd%dat(1:3,j,i)
	  end do
	 end do
	close(1)

return
end

end module
