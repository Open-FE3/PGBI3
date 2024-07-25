# Coded on Oct. 27, 2015
vhfile="./vxyz/vh"${1}".dat"
vzfile="./vxyz/vz"${1}".dat"
out21="./vxyz/vxyz"${1}".ps"
minute=`echo "scale=2; ${i}/6" | bc | sed -e 's/^\./0./g' `
pos5file="../gebco/pos5.dat"
topofile="../gebco/topo.xyz"
polygonfile="../gebco/polygon2.dat"
polygongmt="polygon2.gmt"
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
program tmp
implicit none
integer(4) :: npoly,i,j,k,l,m,n,iclose
real(8) :: cx1,cy1,cx,cy
!###
integer(4) :: nodes,nodeki,node_ki, node,nele,surfptr
real(8),allocatable,dimension(:,:) :: xyz
real(8),allocatable,dimension(:) :: bz
!###################  file read start #####################
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
gfortran tmp.f90
./a.out
########################   get land bz
grdvz="vz.grd"
topogrd="topo.grd"
CPT="vz.cpt"
WESN=-1500/1500/-1500/1500
scl=15/15
bb=a300:"distance\(km\)":/a300:"distance\(km\)":WeSn
makecpt -Cpolar -T-0.1/0.1/0.01 -V > ${CPT}
surface ${vzfile} -G$grdvz -I2/2 -T1 -R$WESN
surface $topofile -G$topogrd -I2/2 -T1 -R$WESN
grdimage $grdvz -B$bb -C${CPT} -JX${scl} -K -R${WESN} -V  -X3 -Y4 > $out21
#grdcontour $topogrd -C1000  -L-7100/-6900 -W0.1/0/255/0  -JX -V -K -O >> $out21
psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -V -W2/0/0/0 >> $out21
psscale -D13/5/8/0.3 -B3  -C$CPT -O -V >> $out21
rm $CPT
rm $grdf $topogrd
rm $polygongmt
gv $out21 &
