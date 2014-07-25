datafile = "statistics.dat"
prefix = ""
number = system("ls positions*png 2> /dev/null | wc -l | bc")

set samples 5000

set terminal pngcairo size 2000,2000 enhanced color font "CMU Sans Serif,28" lw 3
set output prefix."density".number.".png"

set xlabel "z Coordinate"
set ylabel "number density"
plot datafile u 1:2 i 0 smooth csplines title "3 patches" lc 3,\
     datafile u 1:2:(0.2) i 0 w circles notitle lc 3 fs solid,\
     datafile u 1:3 i 0 smooth csplines title "2 patches" lc 1,\
     datafile u 1:3:(0.2) i 0 w circles notitle lc 1 fs solid

set output prefix."totaldensity".number.".png"

set xlabel "z Coordinate"
set ylabel "total number density"
plot datafile u 1:($2+$3) i 0 smooth csplines notitle lc 3,\
     datafile u 1:($2+$3):(0.2) i 0 w circles notitle lc 3 fs solid

set output prefix."bonds".number.".png"

set xlabel "z Coordinate"
set ylabel "formed bonds / possible bonds"
plot datafile u 1:2 i 1 smooth csplines title "3 patches" lc 3,\
     datafile u 1:2:(0.2) i 1 w circles notitle lc 3 fs solid,\
     datafile u 1:3 i 1 smooth csplines title "2 patches" lc 1,\
     datafile u 1:3:(0.2) i 1 w circles notitle lc 1 fs solid
