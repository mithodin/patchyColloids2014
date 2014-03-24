#!/bin/bash
num=`ls positions*png 2> /dev/null | wc -l | bc`
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
	open="eom"
elif [[ "$unamestr" == 'Darwin' ]]; then
	open="open"
fi

for datei in `ls positions*dat`
	do
	vars=`echo $datei | cut -d '-' -f 2-`
	mv $datei positions.dat
	mv statistics-$vars statistics.dat
	gnuplot plotStats.gnu
	gnuplot plotPositions.gnu
	mv positions.dat $datei
	mv statistics.dat statistics-$vars
	$open bonds$num.png density$num.png positions$num.png
	num=$(( $num + 1 ))
	echo $num
done

