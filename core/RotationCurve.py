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

    def print_info(self):
        infoutils.print_header("Rotation curve", self.name, self.description, 0)
        infoutils.print_list_summary(1, "radii in arcsec", self.radii(), description=self.description)
        infoutils.print_list_summary(1, "velocities in km/s", self.velocities())

    def plot(self, label):
        # plt.axhline(y=numpy.array(self.velocities()).mean(), color='red', label='mean')
        plt.plot(self.radii(), self.velocities(), '.', label=label)
        plt.errorbar(self.radii(), self.velocities(), yerr=self.delta_velocities(), fmt=None,
                     marker=None, mew=0)
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$V_{" + label + "},\ km/s$")
        plt.legend()


if __name__ == "__main__":
    rc1 = RotationCurve("../data/ngc338/v_stars_ma.dat", "Stellar MA RC",
                       description="Many-many words about such RC\n with author name and observation parameters. \n "
                                   "Maybe some links and refs.")
    rc1.print_info()
