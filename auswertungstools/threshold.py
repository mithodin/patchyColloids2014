#!/usr/bin/env python2.7
from scipy import *
from sys import argv

def pperc(x):
	return (x+2)/(4*x+2)

for datei in argv[1:]:
	data = genfromtxt(datei, comments="#", delimiter="\t")
	lastdiff=data[0,2]-pperc(data[0,3])
	cut = array([0])
	lastbox = [0,0]
	for box in data:
		diff=box[2]-pperc(box[3])
		if(diff*lastdiff<0):
			cut=append(cut,0.5*(box[1]+lastbox[1]))
		lastdiff=diff
		lastbox=box
	print(str(cut))
	print(datei+": "+str(max(cut)))
