__author__ = 'amarch'
# -*- coding: utf-8 -*-

from bisect import bisect
import sys
import numpy
from scipy import interpolate
from core.velocityEllipsoidReconstr import *
from utils import strutils as infoutils
from utils.bcolors import *
import matplotlib.pyplot as plt

class VelocityDispersion():

    def __init__(self, path, name, description=None):
        self.name = name
        self.path = path
        self.description = description
        self.data_points = []
        self.fake_data_points = []
        self.poly_fit = poly1d([0])
        if os.path.isfile(path):
            sig_file = open(path)
            for line in sig_file:
                if line[0] == '#':
                    pass
                else:
                    line = filter(lambda x: x != '', line.split("  "))
                    if line.__len__() == 3:
                        self.data_points.append((float(line[0]), float(line[1]), float(line[2])))
                    else:
                        self.data_points.append((float(line[0]), float(line[3]), float(line[4])))
            sig_file.close()

    def radii(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in dispersion curve")
        else:
            return map(lambda x: x[0], self.data_points)

    def dispersions(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in dispersion curve")
        else:
            return map(lambda x: x[1], self.data_points)

    def delta_dispersions(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in dispersion curve")
        else:
            return map(lambda x: x[2], self.data_points)

    def add_fake_exp_points(self, point, count, expscale):
        fake_radii = arange(point[0], point[0]+count, 1)
        self.fake_data_points = self.fake_data_points + zip(fake_radii, [point[1]*math.exp(-(x-point[0]) / expscale) for x in fake_radii], [point[2]]*count)

    def expand_minor(self, incl):
        '''Expand minor Sigma(Rcos(i)) -> Sigma(R).
        incl in deg'''
        expand_fun = lambda x: (x[0]/math.cos(incl*math.pi/180), x[1], x[2])
        self.data_points = map(expand_fun, self.data_points)

    def print_info(self, indent):
        infoutils.print_header("Velocity Dispersion", self.name, self.description, indent)
        infoutils.print_simple_param(indent+1, "fake data points count",  str(self.fake_data_points.__len__()))
        if self.poly_fit != poly1d([0]):
            infoutils.print_simple_param(indent+1, "polyfit degree",  self.poly_fit.coeffs.__len__()-1)
            infoutils.print_simple_param(indent+1, "polyfit coeffs",  list(self.poly_fit.coeffs))
        infoutils.print_list_summary(1+indent, "radii in arcsec", self.radii(), description=self.description)
        infoutils.print_list_summary(1+indent, "dispersions in km/s", self.dispersions())


    def plot(self, label, color='red'):
        plt.plot(self.radii(), self.dispersions(), '.', label=label, color=color)
        plt.errorbar(self.radii(), self.dispersions(), yerr=self.delta_dispersions(), fmt='.',
                     marker='.', mew=0, color=color)
        plt.xlabel("$R,\ arcsec$")
        plt.ylabel("$\sigma,\ km/s$")
        if self.poly_fit != poly1d([0]):
            max_R = max(map(abs, self.radii()))
            poly_radii = list(arange(-0.1, -max_R, -0.1)[::-1]) + \
                         list(arange(0, max_R, 0.1))
            poly_sig = list(self.poly_fit(arange(0.1, max_R, 0.1))[::-1]) + list(self.poly_fit(arange(0, max_R, 0.1)))
            plt.plot(poly_radii, poly_sig, '-', label='polinom approx for '+ label)
            if max(poly_sig) > 250:
                plt.ylim(0, 250)
            elif min(poly_sig) < 0:
                plt.ylim(0)
        plt.plot(zip(*(self.fake_data_points))[0], zip(*(self.fake_data_points))[1], 'o', label = "fake points")
        plt.legend(loc='upper right').draw_frame(False)


    def plot_on_one_side(self, label, color='red'):
        plt.plot(map(abs, self.radii()), self.dispersions(), '.', label=label, color=color)
        plt.errorbar(map(abs,self.radii()), self.dispersions(), yerr=self.delta_dispersions(), fmt='.',
                     marker='.', mew=0, color=color)
        plt.xlabel("$R,\ arcsec$")
        plt.ylabel("$\sigma,\ km/s$")