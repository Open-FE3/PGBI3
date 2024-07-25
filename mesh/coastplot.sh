#!/bin/bash
infile="coast.dat"
outfile="coast.ps"
gnuplot <<EOF
set terminal postscript enhanced
set output "$outfile"
set size 0.8,1
set yrange [-2000:2000]
set xrange [-2000:2000]
set grid xtics
set grid ytics 
plot "$infile" u 2:1 w p pt 13 ps 0.2
EOF
gv $outfile &