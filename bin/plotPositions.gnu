coldiameter = 1.0
patchdiameter = 0.11965683746373795115
#patchdiameter = 0.2
datafile = "positions.dat"
params = system("head -n 1 ".datafile." | cut -d '#' -f 2 | tr '\t' '\n'");
height = system("head -n 1 ".datafile." | cut -s -f 7 | cut -d '=' -f 2");
width = system("head -n 1 ".datafile." | cut -s -f 8 | cut -d '=' -f 2");
prefix = ""

set samples 5000
set xrange [-5:5+width]
set yrange [-5:5+height]
set size ratio height/width

set key outside vertical top left reverse
set xlabel "x Coordinate"
set ylabel "z Coordinate"

number = system("ls positions*png 2> /dev/null | wc -l | bc")
set terminal pngcairo size 3000,4000 enhanced color font "CMU Sans Serif,28" lw 3
set output prefix."positions".number.".png"

set object 1 rect from 0,0 to width,height fillstyle empty border lc rgb "gray" lw 2.5

plot datafile u 1:2:(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+cos($3)/2):($2+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+$4*cos($3+2.0/3.0*pi)/2):($2+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+$4*cos($3+4.0/3.0*pi)/2):($2+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+(1-$4)*cos($3+pi)/2):($2+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder\
     	title params
