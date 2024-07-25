#!/bin/sh
#
#
for j in `seq 1 4`
do
 for i in `seq -w 1 1800`
 do
 echo "i=" ${i} "j=" ${j}
  if [ ${j} -le 3 ]; then
   ./bxyzplot_xy.sh 00${i} ${j} # 1 for bx, 2 for by, 3 for bz
  else
   ./vxplot.sh 00${i}
  fi
 done
 if [ ${j} -le 3 ]; then
  convert -delay 30 -loop 1 ./bxyz/bxyz_xy2D00[01][0-9][0-9][0-9].ps bxyz${j}.gif
 else
  convert -delay 30 -loop 1 ./vxyz/vx000[0-4][0-9][0-9].ps vx.gif
 fi
done