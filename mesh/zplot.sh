#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
infile="polygonz.dat"
#infile="gebcoz.dat"
out21="polygonz.ps"
WESN=-500/500/-500/500
scl=15/15
bb=a100:"distance\(km\)":/a100:"distance\(km\)":WeSn
grdf="polygonz.grd"
CPT="bathy.cpt"
#makecpt -Crainbow -T-20/20/0.1 -V > ${CPT}
makecpt -Crainbow -T-9/1/0.1 -V > ${CPT}
#xyz2grd Bathymetry.xyz -GBathymetry.grd -R0/5.488/0/3.402 -I0.014
surface $infile -G$grdf -I2/2 -T1 -R$WESN
grdimage $grdf -B$bb -C${CPT} -JX${scl} -K -R${WESN} -V > $out21
psscale -D17/-1/8/0.3h -B1 -C$CPT -O -V >> $out21
rm $CPT
rm $grdf
gv $out21 &
