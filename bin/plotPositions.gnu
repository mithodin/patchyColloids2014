height = 50
width = 40
coldiameter = 1.0
patchdiameter = 0.119
datafile = "positions-T10.dat"
params = system("head -n 1 ".datafile." | cut -d '#' -f 2 | tr '\t' '\n'");
prefix = "T10"

set samples 5000
set xrange [-1:1+width]
set yrange [-1:1+height]
set size ratio 1

set key outside vertical top left reverse
set xlabel "x Coordinate"
set ylabel "z Coordinate"

number = system("ls positions*png 2> /dev/null | wc -l | bc")
set terminal pngcairo size 2400,2000 enhanced color font "CMU Sans Serif,28"
set output prefix."positions".number.".png"

set object 1 rect from 0,0 to width,height fillstyle empty border lc rgb "gray" lw 2

plot datafile u 1:2:(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+cos($3)/2):($2+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+$4*cos($3+2.0/3.0*pi)/2):($2+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+$4*cos($3+4.0/3.0*pi)/2):($2+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+(1-$4)*cos($3+pi)/2):($2+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width):2:(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+cos($3)/2):($2+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+$4*cos($3+2.0/3.0*pi)/2):($2+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+$4*cos($3+4.0/3.0*pi)/2):($2+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+(1-$4)*cos($3+pi)/2):($2+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width):2:(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+cos($3)/2):($2+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+$4*cos($3+2.0/3.0*pi)/2):($2+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+$4*cos($3+4.0/3.0*pi)/2):($2+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+(1-$4)*cos($3+pi)/2):($2+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width):($2-height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+cos($3)/2):($2-height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+$4*cos($3+2.0/3.0*pi)/2):($2-height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+$4*cos($3+4.0/3.0*pi)/2):($2-height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+(1-$4)*cos($3+pi)/2):($2-height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width):($2-height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+cos($3)/2):($2-height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+$4*cos($3+2.0/3.0*pi)/2):($2-height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+$4*cos($3+4.0/3.0*pi)/2):($2-height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+(1-$4)*cos($3+pi)/2):($2-height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width):($2+height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+cos($3)/2):($2+height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1-width+$4*cos($3+2.0/3.0*pi)/2):($2+height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+$4*cos($3+4.0/3.0*pi)/2):($2+height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1-width+(1-$4)*cos($3+pi)/2):($2+height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width):($2+height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+cos($3)/2):($2+height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+width+$4*cos($3+2.0/3.0*pi)/2):($2+height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+$4*cos($3+4.0/3.0*pi)/2):($2+height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+width+(1-$4)*cos($3+pi)/2):($2+height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u 1:($2+height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+cos($3)/2):($2+height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+$4*cos($3+2.0/3.0*pi)/2):($2+height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+$4*cos($3+4.0/3.0*pi)/2):($2+height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+(1-$4)*cos($3+pi)/2):($2+height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u 1:($2-height):(coldiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+cos($3)/2):($2-height+sin($3)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle, \
     datafile u ($1+$4*cos($3+2.0/3.0*pi)/2):($2-height+$4*sin($3+2.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+$4*cos($3+4.0/3.0*pi)/2):($2-height+$4*sin($3+4.0/3.0*pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder notitle,\
     datafile u ($1+(1-$4)*cos($3+pi)/2):($2-height+(1-$4)*sin($3+pi)/2):(patchdiameter/2.0):($4*2+1) with circles lc variable fs solid noborder\
     	title params
