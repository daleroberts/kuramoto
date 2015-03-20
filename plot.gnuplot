#!/usr/bin/env gnuplot

reset

set terminal pngcairo size 640,480 enhanced font 'Verdana,9'
set output outfile

set style fill transparent solid 0.2 noborder
set border linewidth 1.5

# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75

# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

set title
set xlabel "Time"
set ylabel "Order parameter"

set yrange [0:1.15]

set key font ",8"
set key top left Left reverse samplen -1 

plot infile using 1:($2+$3):($2-$3) with filledcu lc rgb "forest-green" notitle, infile using 1:2 with lines lc rgb "forest-green" title columnheader(1)
