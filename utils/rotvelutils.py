__author__ = 'amarch'
# -*- coding: utf-8 -*-

import numpy as np

def incline_velocity(v, angle):
    return v / np.sin(angle * np.pi / 180.)

def correct_rotation_curve(rdata, vdata, dvdata, r0, v0, incl):
	''' Переносит центр в (r0,v0) и перегибает кривую вращения, 
	а также исправляет за наклон если необходимо.
	'''
	rdata_tmp = [np.abs(r-r0) for r in rdata]
	vdata_tmp = [incline_velocity(abs(v-v0), incl) for v in vdata]
	data = zip(rdata_tmp, vdata_tmp, dvdata)
	data.sort()
	return zip(*data)