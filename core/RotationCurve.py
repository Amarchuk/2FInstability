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


class RotationCurve():

    def __init__(self, path, name, description=None):
        self.name = name
        self.path = path
        self.description = description
        self.data_points = []
        self.fake_data_points = []
        self.poly_fit = poly1d([0])
        if os.path.isfile(path):
            rc_file = open(path)
            for line in rc_file:
                if line[0] == '#':
                    pass
                else:
                    line = filter(lambda x: x != '', line.split("  "))
                    try:
                        self.data_points.append((float(line[0]), float(line[1]), float(line[2])))
                    except ValueError:
                        self.data_points.append((float(line[0]), float(line[1]), 0))
            rc_file.close()

    def radii(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[0], self.data_points)

    def velocities(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[1], self.data_points)

    def delta_velocities(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[2], self.data_points)

    def add_fake_points(self, point, count):
        self.fake_data_points = self.fake_data_points + zip(arange(point[0], point[0]+count, 1), [point[1]]*count, [point[2]]*count)

    def print_info(self, indent):
        infoutils.print_header("Rotation curve", self.name, self.description, indent)
        infoutils.print_simple_param(indent+1, "fake data points count",  str(self.fake_data_points.__len__()))
        if self.poly_fit != poly1d([0]):
            infoutils.print_simple_param(indent+1, "polyfit degree",  self.poly_fit.coeffs.__len__())
            infoutils.print_simple_param(indent+1, "polyfit coeffs",  list(self.poly_fit.coeffs))
        infoutils.print_list_summary(1+indent, "radii in arcsec", self.radii(), description=self.description)
        infoutils.print_list_summary(1+indent, "velocities in km/s", self.velocities())


    def plot(self, label, color='red'):
        plt.plot(self.radii(), self.velocities(), '.', label=label, color=color)
        plt.errorbar(self.radii(), self.velocities(), yerr=self.delta_velocities(), fmt='.',
                     marker='.', mew=0, color=color)
        plt.xlabel("$R,\ arcsec$")
        plt.ylabel("$V,\ km/s$")
        if self.poly_fit != poly1d([0]):
            poly_radii = arange(0, max(self.radii()) + 30, 0.1)
            poly_vel = self.poly_fit(poly_radii)
            plt.plot(poly_radii, poly_vel, '-', label='polinom approx for '+ label)
            if max(poly_vel) > 500:
                plt.ylim(0, 600)
            elif min(poly_vel) < 0:
                plt.ylim(0)
        plt.legend(loc='lower right').draw_frame(False)