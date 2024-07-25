# Coded on May 18, 2016
# This script make a movie composed of surface elevation, uh, and bz
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
width=1500
WESN=-${width}/${width}/-${width}/${width}
bxyzfile="./exyz/exyz_sbtm"${1}".dat"
out21="./exyz/exyz_sbtm"${1}".ps"
pos5file="../gebco/pos5.dat"
topofile="../gebco/topo.xyz"
polygonfile="../gebco/polygon2.dat"
polygongmt="polygon2.gmt"
# get dt
dt=`grep "dt =" ./timepara.f | awk -F[=\.] '{printf($2)}'`
echo "dt=" $dt "[sec]"
minute=` echo "scale=1; ${1} * 10 / 60 " | bc `
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
program tmp
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: ntrik
real(8),allocatable,dimension(:,:) :: xy,exyz,z

open(1,file="./exyz/coord_sbtm.dat")
read(1,*) ntrik
allocate(xy(2,ntrik))
do i=1,ntrik
 read(1,*)xy(1,i),xy(2,i),z
end do
close(1)

allocate(exyz(3,ntrik))
open(1,file="${bxyzfile}")
do i=1,ntrik
 read(1,*) exyz(1,i),exyz(2,i),exyz(3,i)
end do
close(1)

open(1,file="tmp.dat")
do i=1,ntrik
 write(1,*) xy(1,i),xy(2,i),exyz(3,i)
end do
close(1)
!########################################  calculate coastline file
open(1,file="$polygonfile")
open(2,file="$polygongmt")
l=0
100 continue
read(1,*,end=99) l,npoly,iclose
read(1,*) j,cx1,cy1
write(2,'(a4)') "> -Z"
write(2,*) cy1,cx1
do i=2,npoly
read(1,*) j,cx,cy
write(2,*) cy,cx
end do
if (iclose .eq. 1)  write(2,*) cy1,cx1
goto 100
99  continue
if ( l .eq. 0 )then
write(*,*) "GEGEGE there are no lines in the input file ","$infile"
stop
end if
close(1)
close(2)

end program tmp
EOF
########################    tmp.f end  #####################
gfortran tmp.f90
./a.out ## make "tmp.yzc"
#============ bzplot ===================
scl=15/15 ; range=0/15/0/15
bb=a400/a400:"distance\(km\)":WeSn
grdf="bz.grd"
grdtopo="topo.grd"
CPT="bz.cpt"
makecpt -Cpolar -T-5/5/0.05 -V > ${CPT}
surface "tmp.dat" -G$grdf -I2/2 -T1 -R$WESN
surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN
grdimage $grdf -B$bb -C${CPT} -JX${scl} -K -R${WESN} -V  -X3 -Y3 > $out21
psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -V -W2/0/0/0 >> $out21
psxy "$pos5file" -JX$scl -R$WESN -m -K -O -V -L -W2/255/0/0 >> $out21
grdcontour $grdtopo -C1000 -L-8000/-7000 -W0.2/0/255/0 -JX -K -O -V >> $out21
psscale -D15.5/6/10/0.3 -B1 -C$CPT -K -O -V >> $out21
scl2=17/15 ; range2=0/20/0/15
pstext -JX$scl -R$range2  -W0 -G255 -O -V <<EOF >> $out21
16.8   14.5 16 0 4 RM t = $minute [min]
16.8   13.8    16 0 4 RM    ${1} * 10 [s]
17.8 10.3 25 0 4 CM ey
17.8  9.5 16 0 4 CM [nT]
EOF
#=========  plot end =======================
#rm ocaen.dat
rm $CPT a.out nwplocation.dat
rm $polygongmt $grdtopo
rm vec.dat uh.grd uhxyz.dat uh.cpt eta.grd eta.cpt
rm $grdf  #uhxyz.dat
gv $out21 &
