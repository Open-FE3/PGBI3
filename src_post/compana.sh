#!/bin/bash

infile="./bxyz/ORG0_bxyz_ts.dat"
anafile="ana_bxyz.dat"
outfile="./comp_bxyz.ps"

gnuplot <<EOF
set terminal postscript enhanced color
set output "$outfile"
set xlabel "time [min]"
set grid ytics
set grid xtics
set xrange [0:30]
plot "$infile" u (\$1/60):2 w l lw 2 lc rgb "green"  title "cal-bx",\
     "$infile" u (\$1/60):4 w l lw 2 lc rgb "purple" title "cal-bz",\
     "$anafile" u (\$1/60):2 w l lw 2 lc rgb "red"    title "ana-bx",\
     "$anafile" u (\$1/60):4 w l lw 2 lc rgb "blue"    title "ana-bz"
EOF
#plot "$infile" u (\$1/60):2 w p pt 4 lc rgb "green"  title "cal-bx",\
#     "$infile" u (\$1/60):4 w p pt 4 lc rgb "purple" title "cal-bz",\

gv $outfile
