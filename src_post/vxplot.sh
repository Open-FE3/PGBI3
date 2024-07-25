# Coded on June 23, 2016
#
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
width=150
WESN=-${width}/${width}/-${width}/${width}
vxfile="./vxyz/vx"${1}".dat"
out21="./vxyz/vx"${1}".ps"
# get dt
dt=`grep "dt =" ./timepara.f | awk -F[=\.] '{printf($2)}'`
echo "dt=" $dt "[sec]"
minute=` echo "scale=1; ${1} * 10 / 60 " | bc `
#============ vxplot ===================
scl=15/15 ; range=0/15/0/15
bb=a100/a100:"distance\(km\)":WeSn
grdf="bz.grd"
grdtopo="topo.grd"
CPT="bz.cpt"
makecpt -Cpolar -T-40/40/0.1 -V > ${CPT}
surface "$vxfile" -G$grdf -I2/2 -T1 -V -R$WESN
#surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN
grdimage $grdf -B$bb -C${CPT} -JX${scl} -K -R${WESN} -V  -X3 -Y3 > $out21
#psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -V -W2/0/0/0 >> $out21
#psxy "$pos5file" -JX$scl -R$WESN -m -K -O -V -L -W2/255/0/0 >> $out21
#grdcontour $grdtopo -C1000 -L-8000/-7000 -W0.2/0/255/0 -JX -K -O -V >> $out21
psscale -D15.5/6/10/0.3 -B10 -C$CPT -K -O -V >> $out21
psxy  -R -J -K -O -W0 <<EOF >> $out21
-$width 0
$width 0
EOF
scl2=17/15 ; range2=0/20/0/15
pstext -JX$scl2 -R$range2  -W0 -G255 -O -V <<EOF >> $out21
16.8   14.5 16 0 4 RM t = $minute [min]
16.8   13.8    16 0 4 RM    ${1} * 10 [s]
19 13 25 0 4 CM vx
19  12.2 16 0 4 CM [mm/s]
EOF
#=========  plot end =======================
#rm ocaen.dat
rm $CPT a.out
#rm $polygongmt $grdtopo
#rm vec.dat uh.grd uhxyz.dat uh.cpt eta.grd eta.cpt
rm $grdf  tmp.dat
#gv $out21 &
