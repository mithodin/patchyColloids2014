#!/usr/bin/env python2.7
from scipy import *
from sys import argv,stdout
from networkx import *

sigma = 1.0
delta = 0.11965683746373795115
maxd = int(ceil((100**2 + 200**2)**0.5))
bonds = zeros(maxd)
n = zeros(maxd)

gra = empty_graph()

class Particle:
	def __init__(self,x,z,a,kind,name):
		self.x = x
		self.z = z
		self.a = a
		self.k = kind
		self.n = name
	def patchPos(self,i):
		if int(self.k) == 0:
			i = (i%2)*3
		return {
			0: (self.x+cos(self.a)/2.0, self.z+sin(self.a)/2.0),
			1: (self.x+cos(self.a+2.0/3.0*pi)/2.0, self.z+sin(self.a+2.0/3.0*pi)/2.0),
			2: (self.x+cos(self.a+4.0/3.0*pi)/2.0, self.z+sin(self.a+4.0/3.0*pi)/2.0),
			3: (self.x+cos(self.a+pi)/2.0, self.z+sin(self.a+pi)/2.0)
		}[i]
	def bonded(self,p):
		d = ((self.x - p.x)**2 + (self.z-p.z)**2)**0.5
		if d < sigma+delta and d > sigma:
			for i in range(2+int(self.k)):
				for j in range(2+int(p.k)):
					(x1,y1) = self.patchPos(i)
					(x2,y2) = p.patchPos(j)
					if ((x1-x2)**2 + (y1-y2)**2)**0.5 < delta:
						return True
		return False

for datei in argv[1:]:
	data = genfromtxt(datei,comments="#",delimiter="\t")
	colloids = array([])
	for idx, partikel in enumerate(data):
		p = Particle(partikel[0],partikel[1],partikel[2],partikel[3],idx)
		gra.add_node(p.n)
		for p2 in colloids:
			if p.bonded(p2):
				gra.add_edge(p.n,p2.n)
				stdout.write('.')
				stdout.flush()
		colloids = append(colloids,p)
	print ""
	print "All particles added. Scanning..."
	for idx,p1 in enumerate(colloids):
		for p2 in colloids[(idx+1):]:
			d = int(floor(((p1.x-p2.x)**2 + (p1.z-p2.z)**2)**0.5))
			n[d] += 1
			if has_path(gra,p1.n,p2.n):
				bonds[d] += 1
		stdout.write('.')
		stdout.flush()
	for i in range(maxd):
		if n[i] != 0:
			bonds[i] /= n[i]
	savetxt(datei.split(".")[0]+".csv",bonds)
