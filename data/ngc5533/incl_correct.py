#!/usr/bin/python

import sys
import math

incl = 53
radian = incl * math.pi / 180.0
inp = open('s_stars_maN.dat','r')
r_sig = []
sig = []
dsig = []
for line in inp:
	if line[0] == '#':
		print line[:-1]
	else:
		line = filter(lambda x: x != '', line[:-1].split(" "))
		r_sig.append(float(line[0]))
		sig.append(float(line[1]))
		dsig.append(float(line[2]))
		#print float(line[0]), (float(line[1])*math.sin(radian)), (float(line[2])*math.sin(radian)), float(line[3]), float(line[4])
		
inp.close()

inp = open('v_stars_noord.dat','r')
r = []
v = []
dv = []
for line in inp:
	if line[0] == '#':
		print line[:-1]
	else:
		line = filter(lambda x: x != '', line[:-1].split(" "))
		r.append(float(line[0]))
		v.append((float(line[1])*math.sin(radian)))
		dv.append(float(line[2])*math.sin(radian))
		#print float(line[0]), (float(line[1])*math.sin(radian)), (float(line[2])*math.sin(radian)), float(line[3]), float(line[4])
		
inp.close()

for d in r[::-1]:
	ind1 = r.index(d)
	if r_sig.count(d) > 0:
		ind2 = r_sig.index(d)
		print -1*d,-v[ind1],dv[ind1],sig[ind2],dsig[ind2]
	else:
		less = filter(lambda x: x < d, r_sig)
		gt = filter(lambda x: x > d, r_sig)
		if gt.__len__() == 0:
			d1 = less[-2]
			d2 = less[-1]
		else:
			d1 = less[-1]
			d2 = gt[0]
		ind2 = r_sig.index(d1)
		ind3 = r_sig.index(d2)
		s = sig[ind2] + (sig[ind3]-sig[ind2])*(d-d1)/(d2-d1)
		ds = dsig[ind2] + (dsig[ind3]-dsig[ind2])*(d-d1)/(d2-d1)
		print (-1*d),-v[ind1],dv[ind1],round(s,1),round(ds,1)


for d in r:
	ind1 = r.index(d)
	if r_sig.count(d) > 0:
		ind2 = r_sig.index(d)
		print d,v[ind1],dv[ind1],sig[ind2],dsig[ind2]
	else:
		less = filter(lambda x: x < d, r_sig)
		gt = filter(lambda x: x > d, r_sig)
		if gt.__len__() == 0:
			d1 = less[-2]
			d2 = less[-1]
		else:
			d1 = less[-1]
			d2 = gt[0]
		ind2 = r_sig.index(d1)
		ind3 = r_sig.index(d2)
		s = sig[ind2] + (sig[ind3]-sig[ind2])*(d-d1)/(d2-d1)
		ds = dsig[ind2] + (dsig[ind3]-dsig[ind2])*(d-d1)/(d2-d1)
		print d,v[ind1],dv[ind1],round(s,1),round(ds,1)


