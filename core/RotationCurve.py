__author__ = 'amarch'
# -*- coding: utf-8 -*-

from bisect import bisect
import sys
import numpy
from scipy import interpolate
from core.velocityEllipsoidReconstr import *
from utils import strutils as infoutils
from utils.bcolors import *

class RotationCurve():

    name = None
    path = "."
    description = None
    data_points = []

    def __init__(self, path, name, description = None):
        self.name = name
        self.path = path
        self.description = description
        rc_file = open(path)
        for line in rc_file:
            if line[0] == '#':
                pass
            else:
                line = filter(lambda x: x != '', line.split("  "))
                self.data_points.append((float(line[0]),float(line[1]),float(line[2])))
        rc_file.close()

    def getRadiusValues(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[0], self.data_points)

    def getVelocityValues(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[1], self.data_points)

    def getDeltaVelocityValues(self):
        if self.data_points.__len__() == 0:
            raise UnboundLocalError("No points in rotation curve")
        else:
            return map(lambda x: x[2], self.data_points)

    def print_info(self):
        infoutils.print_header("Rotation curve", self.name, self.description, 0)
        infoutils.print_list_summary(1, "radiuses in kpc", self.getRadiusValues(), description=self.description)
        infoutils.print_list_summary(1, "velocities in km/s", self.getVelocityValues())

if __name__ == "__main__":
    rc = RotationCurve("../data/ngc338/v_stars_ma.dat", "Stellar MA RC",
                       description = "Many-many words about such RC\n with author name and observation parameters. \n "
                                     "Maybe some links and refs.")
    rc.print_info()