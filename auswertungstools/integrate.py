#!/usr/bin/env python2.7
from scipy import *
from sys import argv

for datei in argv[1:]:
	data = genfromtxt(datei, comments="#", delimiter="\t")
	total1 = 0
	total2 = 0
	for box in data:
		print box
		total1 += box[1]
		total2 += box[2]
	print total1
	print total2
