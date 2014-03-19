#!/bin/bash
num=`ls positions*png 2> /dev/null | wc -l | bc`
gnuplot plotStats.gnu
gnuplot plotPositions.gnu

open bonds$num.png density$num.png positions$num.png
