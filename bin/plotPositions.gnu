set samples 5000
set xrange [-1:21]
set yrange [-1:21]
set size ratio 1

set xlabel "x Coordinate"
set ylabel "z Coordinate"

number = system("ls positions*png 2> /dev/null | wc -l | bc")
set terminal pngcairo size 800,800 enhanced color
set output "positions".number.".png"

set object 1 rect from 0,0 to 20,20 fillstyle empty

plot "positions.dat" u 1:2:(0.5):($4*2+1) with circles lc variable fs solid noborder notitle, \
     "positions.dat" u ($1+$4*cos($3)/2):($2+$4*sin($3)/2):(0.1):($4*2+1) with circles lc variable fs solid noborder notitle, \
     "positions.dat" u ($1+$4*cos($3+2.0/3.0*pi)/2):($2+$4*sin($3+2.0/3.0*pi)/2):(0.1):($4*2+1) with circles lc variable fs solid noborder notitle,\
     "positions.dat" u ($1+$4*cos($3+4.0/3.0*pi)/2):($2+$4*sin($3+4.0/3.0*pi)/2):(0.1):($4*2+1) with circles lc variable fs solid noborder notitle,\
     "positions.dat" u ($1+(1-$4)*cos($3)/2):($2+(1-$4)*sin($3)/2):(0.1):($4*2+1) with circles lc variable fs solid noborder notitle,\
     "positions.dat" u ($1+(1-$4)*cos($3+pi)/2):($2+(1-$4)*sin($3+pi)/2):(0.1):($4*2+1) with circles lc variable fs solid noborder notitle

unset output
