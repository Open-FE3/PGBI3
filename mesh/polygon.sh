# Coded by T. MINAMI on Jan. 10, 2014
#!/bin/bash
infile=polygon${1}.dat
intermf=polygon${1}.gmt
pos5file="pos5.dat"
topofile="topo.xyz"
out21=polygon${1}.ps
echo $infile
rw=-1300
re=1300
rs=-1300
rn=1300
scale=15/15
bb=a100:"distance\(km\)":/a100:"distance\(km\)":WeSn
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
    program tmp
    implicit none
    integer(4) :: l,npoly,j,i,iclose
    real(8) :: cx1,cy1,cx,cy
!###################  file read start #####################
	open(1,file="$infile")
    open(2,file="$intermf")
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
#############################################
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
WESN=$rw/$re/$rs/$rn
scl=$scale
grdtopo="topo.grd"
#CPT="mesh.cpt"
#makecpt -Cpanoply -T-2.22/-0.2/0.0005 -V > ${CPT}
## note that -L option enclose the polygon automatically in psxy
surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN
psxy "$intermf" -JX$scl -R$WESN  -B"$bb" -K -V -W1,black -X3 -Y4 > $out21
psxy "$pos5file" -JX$scl -R$WESN -O -V -L -W1,black >> $out21 # -L forces "closed"
#grdcontour $grdtopo -C1000 -A1000t -L-7500/-6500 -W0.2 -JX -O -V >> $out21
#pstext "elenum.dat"  -Jx$scl -R$WESN -K -O -V  >> $out21
#pstext  -JX$scl -R$WESN -O -V -W3 <<EOF >> $out21
#2000 0 16 0 4 CM mesh
#EOF
rm tmp.f90
#rm $intermf $grdtopo
#rm mesh.cpt
#rm tmp.yzc
gv $out21 &
