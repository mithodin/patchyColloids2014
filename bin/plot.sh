#!/bin/bash
num=`ls positions*png 2> /dev/null | wc -l | bc`

for datei in `ls positions*dat`
	do
	vars=`echo $datei | cut -d '-' -f 2-3`
	echo $vars
	mv $datei positions.dat
	mv statistics-$vars statistics.dat
	gnuplot plotStats.gnu
	gnuplot plotPositions.gnu
	mv positions.dat $datei
	mv statistics.dat statistics-$vars
	open bonds$num.png density$num.png positions$num.png
	num=$(( $num + 1 ))
	echo $num
done

