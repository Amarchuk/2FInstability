__author__ = 'amarch'
# -*- coding: utf-8 -*-

from utils import strutils as infoutils
import itertools

from scipy.integrate import *

from RotationCurve import *
from Galaxy import *
from utils import strutils as infoutils
import itertools
import copy
from utils import bezier_curve as bezier


class VelocityDispersionHandler():

    '''Line-of-sight velocity dispersion handler.'''

    def __init__(self, galaxyname, galaxy):
        self.name = galaxyname + " Velocity Dispersion Handler"
        self.galaxy = galaxy
        self.sig_ma = None
        self.sig_mi = None

    def set_sig_ma(self, sig_ma):
        self.sig_ma = sig_ma

    def set_sig_mi(self, sig_mi):
        self.sig_mi = sig_mi

    def interpolate_poly_sig(self, sig, deg):
        '''Approximation of folded-over through zero curve (R > 0)'''
        sig.poly_fit = poly1d(polyfit(map(abs, tuple(sig.radii()) + zip(*(sig.fake_data_points))[0]),
                                     tuple(sig.dispersions()) + zip(*(sig.fake_data_points))[1], deg=deg))
                                     # ,w=map(lambda x: 1 / (x + 0.1) ** 2,
                                     #       tuple(sig.delta_dispersions()) + zip(*(sig.fake_data_points))[2])))

    def interpolate_bezier(self, sig, nTimes=1000):
        points = map(lambda p : [abs(p[0]), p[1]], sig.data_points)
        points.sort()
        sig.bezier = bezier.bezier_curve(points, nTimes)


    def print_info(self, indent):
        infoutils.print_header("Velocity Dispersion Handler", self.name, indent=indent)
        self.sig_ma.print_info(indent + 1)
        self.sig_mi.print_info(indent + 1)


    def plot_two_in_one(self):
        plt.title("Line-of-sight dispersions for " + self.name)
        self.sig_ma.plot("$\sigma_{los}^{maj}$")
        self.sig_mi.plot("$\sigma_{los}^{min}$", color="black")


    def plot_two_in_subplots(self):
        f,(ax1,ax2) = plt.subplots(nrows=2, ncols=1)
        plt.subplot(2,1,1)
        plt.title("Line-of-sight dispersions for " + self.name)
        self.sig_ma.plot('$\sigma_{los}^{maj}$')
        plt.subplot(2,1,2)
        self.sig_mi.plot("$\sigma_{los}^{min}$")
        plt.xlim(ax1.get_xlim()[0], ax1.get_xlim()[1])